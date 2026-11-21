% 3-D Cubli 控制对比：PD vs Adaptive vs Fixed-Learning vs FSLC
function reproduce_simulation
%% 1. 系统物理参数初始化
Theta_0 = [16, -5.6, -5.6; -5.6, 16, -5.6; -5.6, -5.6, 16] * 1e-3; 
Theta_w = diag([1.811, 1.811, 1.811]) * 1e-4; 
hat_Theta_0 = Theta_0 - Theta_w; 
m_dist = [0.0787; 0.0787; 0.0787]; 
g_vec = [0; 0; -9.8];              
max_torque = 0.26; % 物理极限力矩             

% 增加一个全局阻尼系数，供所有控制器使用以计算 Bias
K_dump = 0.00005;

%% 2. 目标姿态与仿真设置
v_body = [1; 1; 1] / sqrt(3);
v_world = [0; 0; 1];
axis_rot = cross(v_body, v_world);
axis_rot = axis_rot / norm(axis_rot);
angle_rot = acos(dot(v_body, v_world));
q_target = [cos(angle_rot/2); axis_rot * sin(angle_rot/2)];
dt = 0.001;                        
T_final = 10;                      
time = 0:dt:T_final;
n = length(time);
% 统一的初始姿态
q0 = q_target + [0; 0.03; -0.02; 0.02]; 
q0 = q0 / norm(q0); 
% 有界反作用轮力矩扰动 (四组仿真共用，持续正弦扰动)
dist_all = zeros(3, n);
for i = 1:n
    dist_all(:, i) = [0.02 * sin(2*pi * 1.0 * time(i)); 
                      0.015 * cos(2*pi * 0.5 * time(i)); 
                      0.01 * sin(2*pi * 2.0 * time(i))];
end

%% ================== 仿真 1: FSLC 自学习控制 ==================
x_olc = zeros(10, n); x_olc(4:7, 1) = q0;
u_olc = zeros(3, n);
tilt_err_olc = zeros(2, n); % Roll/Pitch 倾角误差(rad)
Kp_olc = 12.0; Kd_olc = 1.5; 
alpha_olc = 0.85;     
epsilon_olc = 0.05;   
beta_olc = 50000.0;   
sigma_s = 3.0;        
p2_olc = 0.8;         
for i = 1:n-1
    w_h = x_olc(1:3, i); q = x_olc(4:7, i); w_w = x_olc(8:10, i);    
    
    q_inv = [q_target(1); -q_target(2:4)];
    q_err = zeros(4,1);
    q_err(1)   = q_inv(1)*q(1) - dot(q_inv(2:4), q(2:4));
    q_err(2:4) = q_inv(1)*q(2:4) + q(1)*q_inv(2:4) + cross(q_inv(2:4), q(2:4));
    
    Bias = q_err(2:4) + K_dump * w_w; 
    Gyro = w_h;
    tilt_err_olc(:, i) = cubliRollPitchTilt(q, v_body); 
    
    s_surf = Gyro + sigma_s * Bias; 
    s_sq = dot(s_surf, s_surf); 
    s_4 = s_sq * s_sq;          
    
    if i == 1
        u_prev_4 = 0;
    else
        u_prev_sq = dot(u_olc(:, i-1), u_olc(:, i-1));
        u_prev_4 = u_prev_sq * u_prev_sq;
    end
    
    p1_weight = alpha_olc / (1.0 + beta_olc * s_4 + 50.0 * u_prev_4) + epsilon_olc;
    
    u_pd = Kp_olc * Bias + Kd_olc * Gyro; 
    
    if i == 1
        u_raw = p2_olc * u_pd;
    else
        u_raw = p1_weight .* u_olc(:, i-1) + p2_olc * u_pd; 
    end
    u_olc(:, i) = max(min(u_raw, max_torque), -max_torque);
    
    % --- 物理动力学 ---
    p_wh = hat_Theta_0 * w_h + Theta_w * (w_h + w_w);
    R = [1-2*(q(3)^2+q(4)^2), 2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3));
         2*(q(2)*q(3)+q(1)*q(4)), 1-2*(q(2)^2+q(4)^2), 2*(q(3)*q(4)-q(1)*q(2));
         2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)), 1-2*(q(2)^2+q(3)^2)]';
    dot_p_wh = -cross(w_h, p_wh) + cross(m_dist, R * g_vec);
    dot_p_ww = u_olc(:, i) + dist_all(:, i); 
    dot_wh = hat_Theta_0 \ (dot_p_wh - dot_p_ww);
    dot_ww = (Theta_w \ dot_p_ww) - dot_wh;
    W = [0, -w_h(1), -w_h(2), -w_h(3); w_h(1), 0, w_h(3), -w_h(2);
         w_h(2), -w_h(3), 0, w_h(1); w_h(3), w_h(2), -w_h(1), 0];
    dot_q = 0.5 * W * q;
    
    x_olc(1:3, i+1) = x_olc(1:3, i) + dot_wh * dt;
    x_olc(4:7, i+1) = x_olc(4:7, i) + dot_q * dt;
    x_olc(4:7, i+1) = x_olc(4:7, i+1) / norm(x_olc(4:7, i+1)); 
    x_olc(8:10, i+1)= x_olc(8:10, i) + dot_ww * dt;
end
q = x_olc(4:7, n);
tilt_err_olc(:, n) = cubliRollPitchTilt(q, v_body);

%% ================== 仿真 2: 纯自适应控制 (PD 控制) ==================
x_adp = zeros(10, n); x_adp(4:7, 1) = q0;
u_adp = zeros(3, n);
tilt_err_adp = zeros(2, n); % Roll/Pitch 倾角误差(rad)
Kp_base = 8.0; Kd_base = 1.0; Delta_Kp = 15.0; Delta_Kd = 2.0; alpha = 200.0;
for i = 1:n-1
    w_h = x_adp(1:3, i); q = x_adp(4:7, i); w_w = x_adp(8:10, i);    
    
    q_inv = [q_target(1); -q_target(2:4)];
    q_err = zeros(4,1);
    q_err(1)   = q_inv(1)*q(1) - dot(q_inv(2:4), q(2:4));
    q_err(2:4) = q_inv(1)*q(2:4) + q(1)*q_inv(2:4) + cross(q_inv(2:4), q(2:4));
    
    Bias_real = q_err(2:4); 
    Gyro = w_h;
    tilt_err_adp(:, i) = cubliRollPitchTilt(q, v_body); 
    Bias_adp = Bias_real + K_dump * w_w;
    
    e_norm_sq = dot(Bias_adp, Bias_adp); 
    Kp_adapt = Kp_base + Delta_Kp * (1 - exp(-alpha * e_norm_sq));
    Kd_adapt = Kd_base + Delta_Kd * (1 - exp(-alpha * e_norm_sq));
    
    u_raw = Kp_adapt * Bias_adp + Kd_adapt * Gyro; 
    u_adp(:, i) = max(min(u_raw, max_torque), -max_torque);
    
    % --- 物理动力学 ---
    p_wh = hat_Theta_0 * w_h + Theta_w * (w_h + w_w);
    R = [1-2*(q(3)^2+q(4)^2), 2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3));
         2*(q(2)*q(3)+q(1)*q(4)), 1-2*(q(2)^2+q(4)^2), 2*(q(3)*q(4)-q(1)*q(2));
         2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)), 1-2*(q(2)^2+q(3)^2)]';
    dot_p_wh = -cross(w_h, p_wh) + cross(m_dist, R * g_vec);
    dot_p_ww = u_adp(:, i) + dist_all(:, i); 
    dot_wh = hat_Theta_0 \ (dot_p_wh - dot_p_ww);
    dot_ww = (Theta_w \ dot_p_ww) - dot_wh;
    W = [0, -w_h(1), -w_h(2), -w_h(3); w_h(1), 0, w_h(3), -w_h(2);
         w_h(2), -w_h(3), 0, w_h(1); w_h(3), w_h(2), -w_h(1), 0];
    dot_q = 0.5 * W * q;
    
    x_adp(1:3, i+1) = x_adp(1:3, i) + dot_wh * dt;
    x_adp(4:7, i+1) = x_adp(4:7, i) + dot_q * dt;
    x_adp(4:7, i+1) = x_adp(4:7, i+1) / norm(x_adp(4:7, i+1)); 
    x_adp(8:10, i+1)= x_adp(8:10, i) + dot_ww * dt;
end
q = x_adp(4:7, n);
tilt_err_adp(:, n) = cubliRollPitchTilt(q, v_body);

%% ================== 仿真 3: 固定学习强度控制 (PD 控制) ==================
x_fix = zeros(10, n); x_fix(4:7, 1) = q0;
u_fix = zeros(3, n);
tilt_err_fix = zeros(2, n); % Roll/Pitch 倾角误差(rad)
Kp_fix = 12.0; Kd_fix = 1.5; 
p1_fixed = 0.2; % 历史权重
p2_fixed = 0.8; % 当前计算结果权重
for i = 1:n-1
    w_h = x_fix(1:3, i); q = x_fix(4:7, i); w_w = x_fix(8:10, i);    
    
    q_inv = [q_target(1); -q_target(2:4)];
    q_err = zeros(4,1);
    q_err(1)   = q_inv(1)*q(1) - dot(q_inv(2:4), q(2:4));
    q_err(2:4) = q_inv(1)*q(2:4) + q(1)*q_inv(2:4) + cross(q_inv(2:4), q(2:4));
    
    Bias = q_err(2:4) + K_dump * w_w; 
    Gyro = w_h;
    tilt_err_fix(:, i) = cubliRollPitchTilt(q, v_body); 
    
    u_pd = Kp_fix * Bias + Kd_fix * Gyro; 
    
    if i == 1
        u_raw = p2_fixed * u_pd;
    else
        u_raw = p2_fixed * u_pd + p1_fixed * u_fix(:, i-1); 
    end
    u_fix(:, i) = max(min(u_raw, max_torque), -max_torque);
    
    % --- 物理动力学 ---
    p_wh = hat_Theta_0 * w_h + Theta_w * (w_h + w_w);
    R = [1-2*(q(3)^2+q(4)^2), 2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3));
         2*(q(2)*q(3)+q(1)*q(4)), 1-2*(q(2)^2+q(4)^2), 2*(q(3)*q(4)-q(1)*q(2));
         2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)), 1-2*(q(2)^2+q(3)^2)]';
    dot_p_wh = -cross(w_h, p_wh) + cross(m_dist, R * g_vec);
    dot_p_ww = u_fix(:, i) + dist_all(:, i); 
    dot_wh = hat_Theta_0 \ (dot_p_wh - dot_p_ww);
    dot_ww = (Theta_w \ dot_p_ww) - dot_wh;
    W = [0, -w_h(1), -w_h(2), -w_h(3); w_h(1), 0, w_h(3), -w_h(2);
         w_h(2), -w_h(3), 0, w_h(1); w_h(3), w_h(2), -w_h(1), 0];
    dot_q = 0.5 * W * q;
    
    x_fix(1:3, i+1) = x_fix(1:3, i) + dot_wh * dt;
    x_fix(4:7, i+1) = x_fix(4:7, i) + dot_q * dt;
    x_fix(4:7, i+1) = x_fix(4:7, i+1) / norm(x_fix(4:7, i+1)); 
    x_fix(8:10, i+1)= x_fix(8:10, i) + dot_ww * dt;
end
q = x_fix(4:7, n);
tilt_err_fix(:, n) = cubliRollPitchTilt(q, v_body);

%% ================== 仿真 4: PD 控制 (保留硬件调试系数 0.7 和 /10) ==================
x_hw = zeros(10, n); x_hw(4:7, 1) = q0;
u_hw = zeros(3, n);
tilt_err_hw = zeros(2, n); % Roll/Pitch 倾角误差(rad)

% 保留硬件调试系数：比例项乘 0.7，微分项除以 10。
% 保留变量名 x_hw/u_hw 以尽量减少对原绘图结构的修改。
Balance_KP = Kp_olc; 
Balance_KD = Kd_olc;  

for i = 1:n-1
    w_h = x_hw(1:3, i); q = x_hw(4:7, i); w_w = x_hw(8:10, i);    
    
    q_inv = [q_target(1); -q_target(2:4)];
    q_err = zeros(4,1);
    q_err(1)   = q_inv(1)*q(1) - dot(q_inv(2:4), q(2:4));
    q_err(2:4) = q_inv(1)*q(2:4) + q(1)*q_inv(2:4) + cross(q_inv(2:4), q(2:4));
    
    Bias = q_err(2:4) + K_dump * w_w; 
    Gyro = w_h;
    tilt_err_hw(:, i) = cubliRollPitchTilt(q, v_body); 
    
    % PD：按硬件实验控制律计算（以仿真中的 Gyro 单位为准）
    u_raw = Balance_KP * 0.7 * Bias + (Balance_KD / 10.0) * Gyro; 
    
    u_hw(:, i) = max(min(u_raw, max_torque), -max_torque);
    
    % --- 物理动力学 ---
    p_wh = hat_Theta_0 * w_h + Theta_w * (w_h + w_w);
    R = [1-2*(q(3)^2+q(4)^2), 2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3));
         2*(q(2)*q(3)+q(1)*q(4)), 1-2*(q(2)^2+q(4)^2), 2*(q(3)*q(4)-q(1)*q(2));
         2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)), 1-2*(q(2)^2+q(3)^2)]';
    dot_p_wh = -cross(w_h, p_wh) + cross(m_dist, R * g_vec);
    dot_p_ww = u_hw(:, i) + dist_all(:, i); 
    dot_wh = hat_Theta_0 \ (dot_p_wh - dot_p_ww);
    dot_ww = (Theta_w \ dot_p_ww) - dot_wh;
    W = [0, -w_h(1), -w_h(2), -w_h(3); w_h(1), 0, w_h(3), -w_h(2);
         w_h(2), -w_h(3), 0, w_h(1); w_h(3), w_h(2), -w_h(1), 0];
    dot_q = 0.5 * W * q;
    
    x_hw(1:3, i+1) = x_hw(1:3, i) + dot_wh * dt;
    x_hw(4:7, i+1) = x_hw(4:7, i) + dot_q * dt;
    x_hw(4:7, i+1) = x_hw(4:7, i+1) / norm(x_hw(4:7, i+1)); 
    x_hw(8:10, i+1)= x_hw(8:10, i) + dot_ww * dt;
end
q = x_hw(4:7, n);
tilt_err_hw(:, n) = cubliRollPitchTilt(q, v_body);

% 0:dt:T_final 共 n 个采样时刻，但控制计算只执行到 n-1。
% 绘图中对末端控制量进行零阶保持，避免 10 s 处人为掉到零。
u_olc(:, n) = u_olc(:, n-1);
u_adp(:, n) = u_adp(:, n-1);
u_fix(:, n) = u_fix(:, n-1);
u_hw(:, n)  = u_hw(:, n-1);

%% ================== PD 与 Fixed-Learning 差异程度检查 ==================
idx_check = time >= 5 & time <= 8;
roll_rmse  = sqrt(mean((tilt_err_hw(1,idx_check) - tilt_err_fix(1,idx_check)).^2));
pitch_rmse = sqrt(mean((tilt_err_hw(2,idx_check) - tilt_err_fix(2,idx_check)).^2));
fprintf('5-8 s PD 与 Fixed-Learning 的 Roll RMSE = %.6e rad\n', roll_rmse);
fprintf('5-8 s PD 与 Fixed-Learning 的 Pitch RMSE = %.6e rad\n', pitch_rmse);

%% ================== 5. 绘图分析：论文双栏版，左右等宽，仅导出 PNG ==================

save_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

figure_names = {'Attitude Error', ...
                'Body Angular Velocity', ...
                'Wheel Angular Velocity', ...
                'Control Torque'};

y_labels_list = {
    {'Roll (rad)', 'Pitch (rad)'}, ...
    {'\omega_{hx} (rad/s)', '\omega_{hy} (rad/s)', '\omega_{hz} (rad/s)'}, ...
    {'\omega_{w1} (rad/s)', '\omega_{w2} (rad/s)', '\omega_{w3} (rad/s)'}, ...
    {'T_1 (N\cdot m)', 'T_2 (N\cdot m)', 'T_3 (N\cdot m)'}
};

% ================== 论文绘图统一参数 ==================
font_name = 'Times New Roman';

axis_font_size   = 12;     % 坐标轴刻度字体
label_font_size  = 15;     % x/y 标签字体
legend_font_size = 10;     % 图例字体

line_width_main  = 1.8;    % 左侧主图线宽
line_width_inset = 1.6;    % 右侧局部放大图线宽
line_width_pd    = 2.3;    % PD 蓝色实线适当加粗，提高可辨识度
axis_line_width  = 1.1;    % 坐标轴线宽

% 使用 MATLAB 默认颜色顺序，但固定各方法的颜色与绘制顺序无关
% PD: 第1色；Adaptive: 第2色；Fixed-Learning: 第3色；FSLC: 第4色
method_colors = get(groot, 'defaultAxesColorOrder');

for fig_idx = 1:4

    % ================== 创建适合双栏论文的画布 ==================
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'centimeters', ...
        'Position', [2, 2, 19.5, 13.8], ...
        'Name', figure_names{fig_idx});

    set(fig, 'ToolBar', 'none', 'MenuBar', 'none');

    if fig_idx == 1
        num_subplots = 2;   % Attitude Error 只有两行
    else
        num_subplots = 3;
    end

    % 不使用总标题
    % sgtitle(figure_names{fig_idx}, 'FontSize', 13, 'FontWeight', 'bold');

    % ================== 左右两列等宽布局 ==================
    % 增大左边距，避免纵坐标标签显示不全
    main_x  = 0.105;
    main_w  = 0.385;

    inset_x = 0.57;
    inset_w = 0.385;

    top_y    = 0.94;
    bottom_y = 0.10;
    row_gap  = 0.075;

    row_h = (top_y - bottom_y - (num_subplots - 1) * row_gap) / num_subplots;

    for j = 1:num_subplots

        % ================== 数据映射 ==================
        % fig_idx=1: 统一从机体对角线相对世界竖直方向计算 Roll/Pitch 倾角误差
        if fig_idx == 1
            data_olc = tilt_err_olc(j,:);
            data_adp = tilt_err_adp(j,:);
            data_fix = tilt_err_fix(j,:);
            data_hw  = tilt_err_hw(j,:);
        elseif fig_idx == 2
            data_olc = x_olc(j,:);
            data_adp = x_adp(j,:);
            data_fix = x_fix(j,:);
            data_hw  = x_hw(j,:);
        elseif fig_idx == 3
            data_olc = x_olc(j+7,:);
            data_adp = x_adp(j+7,:);
            data_fix = x_fix(j+7,:);
            data_hw  = x_hw(j+7,:);
        elseif fig_idx == 4
            data_olc = u_olc(j,:);
            data_adp = u_adp(j,:);
            data_fix = u_fix(j,:);
            data_hw  = u_hw(j,:);
        end

        % 当前行位置
        y_pos = top_y - j * row_h - (j - 1) * row_gap;

        % ==========================================================
        % 左侧主图：0-10 s
        % ==========================================================
        ax = axes(fig, 'Position', [main_x, y_pos, main_w, row_h]);

        % 先画其他三种方法，再把 PD 蓝色实线置顶。
        % 保留所有控制器原始数据，不改变绘图的物理量。
        hold(ax, 'on');
        p2 = plot(ax, time, data_adp, 'Color', method_colors(2,:), ...
                  'LineStyle', '-', 'LineWidth', line_width_main);
        p3 = plot(ax, time, data_fix, 'Color', method_colors(3,:), ...
                  'LineStyle', '-', 'LineWidth', line_width_main);
        p4 = plot(ax, time, data_olc, 'Color', method_colors(4,:), ...
                  'LineStyle', '-', 'LineWidth', line_width_main);
        p1 = plot(ax, time, data_hw, 'Color', [0 0.4470 0.7410], ...
                  'LineStyle', '-', 'LineWidth', line_width_pd);

        grid(ax, 'on');
        box(ax, 'on');
        xlim(ax, [0, 10]);

        ylabel(ax, y_labels_list{fig_idx}{j}, ...
            'FontSize', label_font_size, ...
            'FontWeight', 'bold', ...
            'FontName', font_name, ...
            'Interpreter', 'tex');

        if j == num_subplots
            xlabel(ax, 'Time (s)', ...
                'FontSize', label_font_size, ...
                'FontWeight', 'bold', ...
                'FontName', font_name);
        end

        set(ax, ...
            'FontName', font_name, ...
            'FontSize', axis_font_size, ...
            'LineWidth', axis_line_width);

        % ================== 图例：放在左侧第一幅曲线图内部 ==================
if j == 1
    yl = ylim(ax);
    y_range = yl(2) - yl(1);

    % 为图例在图内留出顶部空白
    if fig_idx == 2
        % Body Angular Velocity 第一行曲线幅值较大，顶部多留空
        ylim(ax, [yl(1), yl(2) + 1.10 * y_range]);
    elseif fig_idx == 4
        ylim(ax, [yl(1), yl(2) + 0.80 * y_range]);
    elseif fig_idx == 3
        ylim(ax, [yl(1), yl(2) + 0.55 * y_range]);
    else
        ylim(ax, [yl(1), yl(2) + 0.45 * y_range]);
    end

    % 图例保持在左侧第一个坐标轴内部，不手动设置 Position
    lgd = legend(ax, [p1, p2, p3, p4], ...
        {'PD', 'Adaptive', 'Fixed-Learning', 'FSLC'}, ...
        'Location', 'northeast', ...
        'FontSize', legend_font_size, ...
        'FontName', font_name, ...
        'Box', 'off', ...
        'Color', 'none');

    % 单列图例，PD 显示为蓝色虚线
    lgd.NumColumns = 1;
end

        % ==========================================================
        % 右侧局部放大图：5-8 s，与左图等宽等高
        % ==========================================================
        ax_in = axes(fig, 'Position', [inset_x, y_pos, inset_w, row_h]);

        % 局部放大图与主图保持颜色、线型和遮挡顺序一致
        hold(ax_in, 'on');
        plot(ax_in, time, data_adp, 'Color', method_colors(2,:), ...
             'LineStyle', '-', 'LineWidth', line_width_inset);
        plot(ax_in, time, data_fix, 'Color', method_colors(3,:), ...
             'LineStyle', '-', 'LineWidth', line_width_inset);
        plot(ax_in, time, data_olc, 'Color', method_colors(4,:), ...
             'LineStyle', '-', 'LineWidth', line_width_inset);
        plot(ax_in, time, data_hw, 'Color', [0 0.4470 0.7410], ...
             'LineStyle', '-', 'LineWidth', line_width_pd);

        xlim(ax_in, [5, 8]);

        idx = (time >= 5) & (time <= 8);

        y_min = min([min(data_hw(idx)), min(data_adp(idx)), ...
                     min(data_fix(idx)), min(data_olc(idx))]);

        y_max = max([max(data_hw(idx)), max(data_adp(idx)), ...
                     max(data_fix(idx)), max(data_olc(idx))]);

        y_margin = max(abs(y_max - y_min) * 0.20, 1e-4);
        ylim(ax_in, [y_min - y_margin, y_max + y_margin]);

        grid(ax_in, 'on');
        box(ax_in, 'on');

        % 不使用右侧小图标题
        % title(ax_in, '5s-8s details', 'FontSize', 8, 'FontWeight', 'normal');

        set(ax_in, ...
            'FontName', font_name, ...
            'FontSize', axis_font_size, ...
            'LineWidth', axis_line_width);

        if j == num_subplots
            xlabel(ax_in, 'Time (s)', ...
                'FontSize', label_font_size, ...
                'FontWeight', 'bold', ...
                'FontName', font_name);
        end

    end

    % ================== 文件名处理 ==================
    file_name = regexprep(figure_names{fig_idx}, '\s+', '_');

    % ================== 仅导出 PNG ==================
    exportgraphics(fig, fullfile(save_dir, [file_name, '.png']), ...
        'Resolution', 600);

    close(fig);
end

methods = ["FSLC"; "Fixed-learning"; "Adaptive"; "PD"];
tilts = {tilt_err_olc, tilt_err_fix, tilt_err_adp, tilt_err_hw};
window = time >= 3 & time <= 10;
stats = zeros(4,4);
for k = 1:4
    angles = tilts{k}(:,window);
    stats(k,:) = [mean(angles(1,:)), var(angles(1,:)), ...
                  mean(angles(2,:)), var(angles(2,:))];
end
result = array2table(stats, 'VariableNames', ...
    {'RollMean_rad','RollVariance_rad2','PitchMean_rad','PitchVariance_rad2'});
result = addvars(result,methods,'Before',1,'NewVariableNames','Method');
writetable(result,fullfile(save_dir,'table1.csv'));
ref = readtable(fullfile(fileparts(mfilename('fullpath')),'reference','table1.csv'));
assert(isequal(string(ref.Method),methods) && all(abs(stats-ref{:,2:5}) < 1e-12,'all'), ...
    'Simulation results differ from the archived values.');
save(fullfile(save_dir,'trajectories.mat'),'time','x_olc','x_fix','x_adp','x_hw', ...
    'u_olc','u_fix','u_adp','u_hw','tilt_err_olc','tilt_err_fix','tilt_err_adp','tilt_err_hw');
end

%% Local function
function roll_pitch = cubliRollPitchTilt(q, v_body)
% q: 单位四元数 [w; x; y; z]，将机体系方向旋转到世界系
% v_body: 机体系中沿平衡角点对角线的单位方向 [1;1;1]/sqrt(3)
% 在期望平衡状态，v_world_axis=[0;0;1]，因此 Roll/Pitch 都为 0。
% 注意：这里计算的是对角线在世界坐标系下的两轴倾角，
% 不是直接把四元数 q_e(2:3) 当作角度，也不是整个机体的欧拉角。
    q = q(:) / norm(q);
    v_body = v_body(:) / norm(v_body);
    qv = q(2:4);
    % 四元数将机体系对角线旋转到世界坐标系
    v_world_axis = v_body + 2*q(1)*cross(qv, v_body) ...
                    + 2*cross(qv, cross(qv, v_body));
    v_world_axis = v_world_axis / norm(v_world_axis);

    % 将对角线相对竖直方向的偏离按世界系 X/Y 方向分解
    roll  = atan2(-v_world_axis(2), v_world_axis(3));
    pitch = atan2(v_world_axis(1), hypot(v_world_axis(2), v_world_axis(3)));
    roll_pitch = [roll; pitch];
end

