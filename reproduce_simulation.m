function reproduce_simulation
%% Plant parameters
Theta_0 = [16, -5.6, -5.6; -5.6, 16, -5.6; -5.6, -5.6, 16] * 1e-3; 
Theta_w = diag([1.811, 1.811, 1.811]) * 1e-4; 
hat_Theta_0 = Theta_0 - Theta_w; 
m_dist = [0.0787; 0.0787; 0.0787]; 
g_vec = [0; 0; -9.8];              
max_torque = 0.26;
K_dump = 0.00005;

%% Balance attitude and disturbance
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
q0 = q_target + [0; 0.03; -0.02; 0.02]; 
q0 = q0 / norm(q0); 
dist_all = zeros(3, n);
for i = 1:n
    dist_all(:, i) = [0.02 * sin(2*pi * 1.0 * time(i)); 
                      0.015 * cos(2*pi * 0.5 * time(i)); 
                      0.01 * sin(2*pi * 2.0 * time(i))];
end

%% Intensity-varying self-learning control (FSLC)
x_olc = zeros(10, n); x_olc(4:7, 1) = q0;
u_olc = zeros(3, n);
tilt_err_olc = zeros(2, n);
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

%% Adaptive control
x_adp = zeros(10, n); x_adp(4:7, 1) = q0;
u_adp = zeros(3, n);
tilt_err_adp = zeros(2, n);
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

%% Fixed-learning control
x_fix = zeros(10, n); x_fix(4:7, 1) = q0;
u_fix = zeros(3, n);
tilt_err_fix = zeros(2, n);
Kp_fix = 12.0; Kd_fix = 1.5; 
p1_fixed = 0.2;
p2_fixed = 0.8;
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

%% PD control
x_hw = zeros(10, n); x_hw(4:7, 1) = q0;
u_hw = zeros(3, n);
tilt_err_hw = zeros(2, n);
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
    
    % The factors 0.7 and 1/10 come from the hardware PD implementation.
    u_raw = Balance_KP * 0.7 * Bias + (Balance_KD / 10.0) * Gyro; 
    u_hw(:, i) = max(min(u_raw, max_torque), -max_torque);
    
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

% Hold the last applied torque at the final plotting time.
u_olc(:, n) = u_olc(:, n-1);
u_adp(:, n) = u_adp(:, n-1);
u_fix(:, n) = u_fix(:, n-1);
u_hw(:, n)  = u_hw(:, n-1);

idx_check = time >= 5 & time <= 8;
roll_rmse  = sqrt(mean((tilt_err_hw(1,idx_check) - tilt_err_fix(1,idx_check)).^2));
pitch_rmse = sqrt(mean((tilt_err_hw(2,idx_check) - tilt_err_fix(2,idx_check)).^2));
fprintf('5-8 s PD 与 Fixed-Learning 的 Roll RMSE = %.6e rad\n', roll_rmse);
fprintf('5-8 s PD 与 Fixed-Learning 的 Pitch RMSE = %.6e rad\n', pitch_rmse);

%% Figures 2-5
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

font_name = 'Times New Roman';
axis_font_size   = 12;
label_font_size  = 15;
legend_font_size = 10;

line_width_main  = 1.8;
line_width_inset = 1.6;
line_width_pd    = 2.3;
axis_line_width  = 1.1;
method_colors = get(groot, 'defaultAxesColorOrder');

for fig_idx = 1:4
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'centimeters', ...
        'Position', [2, 2, 19.5, 13.8], ...
        'Name', figure_names{fig_idx});
    set(fig, 'ToolBar', 'none', 'MenuBar', 'none');

    if fig_idx == 1
        num_subplots = 2;
    else
        num_subplots = 3;
    end

    main_x  = 0.105;
    main_w  = 0.385;
    inset_x = 0.57;
    inset_w = 0.385;
    top_y    = 0.94;
    bottom_y = 0.10;
    row_gap  = 0.075;
    row_h = (top_y - bottom_y - (num_subplots - 1) * row_gap) / num_subplots;

    for j = 1:num_subplots
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

        y_pos = top_y - j * row_h - (j - 1) * row_gap;
        ax = axes(fig, 'Position', [main_x, y_pos, main_w, row_h]);
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

if j == 1
    yl = ylim(ax);
    y_range = yl(2) - yl(1);
    if fig_idx == 2
        ylim(ax, [yl(1), yl(2) + 1.10 * y_range]);
    elseif fig_idx == 4
        ylim(ax, [yl(1), yl(2) + 0.80 * y_range]);
    elseif fig_idx == 3
        ylim(ax, [yl(1), yl(2) + 0.55 * y_range]);
    else
        ylim(ax, [yl(1), yl(2) + 0.45 * y_range]);
    end

    lgd = legend(ax, [p1, p2, p3, p4], ...
        {'PD', 'Adaptive', 'Fixed-Learning', 'FSLC'}, ...
        'Location', 'northeast', ...
        'FontSize', legend_font_size, ...
        'FontName', font_name, ...
        'Box', 'off', ...
        'Color', 'none');
    lgd.NumColumns = 1;
end

        ax_in = axes(fig, 'Position', [inset_x, y_pos, inset_w, row_h]);
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

    file_name = regexprep(figure_names{fig_idx}, '\s+', '_');
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

function roll_pitch = cubliRollPitchTilt(q, v_body)
% The scalar-first quaternion rotates the body diagonal into world coordinates.
% Roll and pitch measure the tilt of that diagonal relative to vertical.
    q = q(:) / norm(q);
    v_body = v_body(:) / norm(v_body);
    qv = q(2:4);
    v_world_axis = v_body + 2*q(1)*cross(qv, v_body) ...
                    + 2*cross(qv, cross(qv, v_body));
    v_world_axis = v_world_axis / norm(v_world_axis);

    roll  = atan2(-v_world_axis(2), v_world_axis(3));
    pitch = atan2(v_world_axis(1), hypot(v_world_axis(2), v_world_axis(3)));
    roll_pitch = [roll; pitch];
end
