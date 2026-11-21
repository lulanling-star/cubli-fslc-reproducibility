# Cubli attitude stabilization with intensity-varying self-learning control

Reproduction materials for *Attitude Stabilization of a Reaction-Wheel Robot via Intensity-Varying Self-Learning Control*. The package contains the nonlinear MATLAB model, the four controllers compared in the manuscript, their simulated trajectories, and the four recorded hardware trials used in Section 5.

| Manuscript result | How to reproduce | Output in results/ |
| --- | --- | --- |
| Fig. 2 and Table 1 | Run reproduce_simulation in MATLAB | Attitude_Error.png and table1.csv |
| Figs. 3–5 | Same MATLAB command | Body_Angular_Velocity.png, Wheel_Angular_Velocity.png and Control_Torque.png |
| Fig. 6 and Table 2 | Run reproduce_hardware.py | measured_tilt.png (replot) and table2.csv |

## Numerical simulation

Open the repository root in MATLAB and run the following command. Tested with base MATLAB R2023a.

~~~matlab
reproduce_simulation
~~~

The function initializes the physical parameters, reference and initial quaternions, disturbance, common torque limit and controller gains. It advances the rigid-body and wheel states by forward Euler at 1 ms, normalizing the quaternion after each step. Table 1 uses 7,001 points over the inclusive interval 3–10 s. Roll and pitch are the tilt components of the body diagonal relative to world vertical. Means are in rad and sample variances are in rad².

In results/trajectories.mat the time vector is in seconds. Columns are time samples. The suffixes olc, fix, adp and hw identify FSLC, fixed-learning, adaptive and simulated PD respectively. For each state array x, rows 1–3 contain body rates (rad/s), 4–7 the scalar-first attitude quaternion and 8–10 the wheel rates (rad/s). The u arrays contain controller torques (N·m), with the last command held for plotting. The tilt_err arrays contain roll and pitch (rad).

The script checks all 16 statistics against reference/table1.csv with an absolute tolerance of 1e-12. One printed entry in the submitted Table 1 differs from the recomputation. The adaptive-controller pitch variance is 7.47939375 × 10⁻⁶ rad² (7.4794 × 10⁻⁶ at five significant figures). The submitted table shows 7.4797 × 10⁻⁶ rad². The method ranking is unchanged.

## Hardware data

From the repository root, run the commands below. The replay was tested with Python 3.11.5 and Matplotlib 3.7.2. It reads the included CSVs directly.

~~~shell
python -m pip install -r requirements.txt
python reproduce_hardware.py
~~~

| Controller | Original DataScope export in data/hardware/ |
| --- | --- |
| PD | pid.csv |
| Adaptive | adaptive.csv |
| Fixed-learning | selflearning.csv |
| FSLC | slc.csv |

DataScope exports have a row per channel and use the first column for the row label. The first row gives one-based sample numbers. CH1 and CH2 are the two measured balancing-channel tilt angles in degrees. For each trial the script selects sample numbers 300–800 inclusive (CSV columns 301–801), matching the plotted 3–8 s interval at 10 ms per sample under the manuscript's indexing convention. Table 2 summarizes fluctuations within each 501-sample record using the arithmetic mean and sample variance (divisor 500). The script checks all 16 values at the four decimal places reported in Table 2.

The embedded/Balance_Control_X.c excerpt documents the deployed X-channel FSLC update described in Section 5.2. It relies on firmware-defined gains. Rebuilding the complete STM32 experiment requires the corresponding board project, controller configurations and calibration settings, which are outside this release.
