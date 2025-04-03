% ELEC 4700A - Filter Response of a Grating with a Gaussian Envelope
% Name: Ayooluwa Owolabi and Student#: 101237575
% Milestone 5: Open_ended_milestone

set(0, 'defaultaxesfontsize', 20)
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesLineWidth', 2)
set(0, 'DefaultFigureWindowStyle', 'docked')

% Define constants
c_c = 299792458;
c_eps_0 = 8.8542149e-12;
c_eps_0_cm = c_eps_0/100;
c_mu_0 = 1 / c_eps_0 / c_c^2;
c_hb = 1.05457266913e-34;
c_h = c_hb * 2 * pi;

% Input parameters (Gaussian pulse)
InputParasL.E0 = 1e5;
InputParasL.we = 0;
InputParasL.t0 = 2e-12;
InputParasL.wg = 2e-13;
InputParasL.phi = 0;
InputParasR = 0;

% Reflection coefficients (no reflections)
RL = 0;
RR = 0;

% Waveguide properties
n_g = 3.5;
vg = c_c / n_g * 1e2;
Lambda = 1550e-9;

plotN = 10;

% Define spatial and time steps
L = 1000e-6 * 1e2;
XL = [0, L];
YL = [-2.5*InputParasL.E0, 2.5*InputParasL.E0];

Nz = 500;
dz = L / (Nz - 1);
dt = dz / vg;
fsync = dt * vg / dz;

Nt = floor(2 * Nz);
tmax = Nt * dt;
t_L = dt * Nz;

z = linspace(0, L, Nz).';

% Define grating properties
kappa0 = 75; % Grating strength
central_freq = L / 2; % Grating center
bandwidth = L / 7; % Controls grating width

% Define the spatially varying grating coupling coefficient κ(z)
kappa = kappa0 * exp(-((z - central_freq) / bandwidth).^2) ... % Gaussian envelope for smooth transition
       .* cos(2 * pi * (z / L) * 10); % Sinusoidal modulation for frequency selectivity

% Limit grating region to the middle third of the waveguide
kappa(z < L/3) = 0;
kappa(z > 2*L/3) = 0;

% Plot κ(z) distribution
figure('Name', 'Grating Profile');
plot(z * 1e4, kappa, 'k', 'LineWidth', 2);
xlabel('z (μm)');
ylabel('κ');
title('Grating Coupling Coefficient κ as a Function of z');
grid on;

% Initialize fields and input/output storage
time = nan(1, Nt);
InputL = nan(1, Nt);
InputR = nan(1, Nt);
OutputL = nan(1, Nt);
OutputR = nan(1, Nt);
Ef = zeros(size(z));
Er = zeros(size(z));

% Source function handles
Ef1 = @SourceFct;
ErN = @SourceFct;

% Initial conditions
t = 0;
time(1) = t;
InputL(1) = Ef1(t, InputParasL);
InputR(1) = ErN(t, InputParasR);
OutputR(1) = Ef(Nz);
OutputL(1) = Er(1);
Ef(1) = InputL(1);
Er(Nz) = InputR(1);

% Gain and detuning parameters (set to zero)
beta_r = 0;
beta_i = 0;
beta = ones(size(z)) * (beta_r + 1i * beta_i);
exp_det = exp(-1i * dz * beta);

% Time marching loop for wave propagation
for i = 2:Nt
    t = dt * (i - 1);
    time(i) = t;
    InputL(i) = Ef1(t, InputParasL);
    InputR(i) = ErN(t, InputParasR);
    
    % Reflection boundary conditions
    Ef(1) = InputL(i) + RL * Er(1);
    Er(Nz) = InputR(i) + RR * Ef(Nz);
    
    % Store previous Ef before updating
    Ef_prev = Ef;
    
    % Wave propagation with grating coupling
    Ef(2:Nz) = fsync * exp_det(1:Nz-1) .* Ef(1:Nz-1) + 1i * dz * kappa(2:Nz) .* Er(2:Nz);
    Er(1:Nz-1) = fsync * exp_det(2:Nz) .* Er(2:Nz) + 1i * dz * kappa(1:Nz-1) .* Ef_prev(1:Nz-1);
    
    % Store output at each time step
    OutputR(i) = Ef(Nz) * (1 - RR);
    OutputL(i) = Er(1) * (1 - RL);
end

% FFT Analysis for frequency-domain representation
fftInputL = fftshift(fft(InputL));
fftOutputR = fftshift(fft(OutputR));
fftOutputL = fftshift(fft(OutputL));
omega = fftshift(wspace(time));
omega_THz = omega / (2 * pi * 1e12);

% FFT Magnitude plot
figure('Name', 'FFT Magnitude Response');
plot(omega_THz, abs(fftInputL), 'm', 'DisplayName', 'Input');
hold on;
plot(omega_THz, abs(fftOutputR), 'k', 'DisplayName', 'Right Output', 'LineWidth', 2);
plot(omega_THz, abs(fftOutputL), 'c', 'DisplayName', 'Left Output', 'LineWidth', 2);
hold off;
xlabel('Frequency (THz)');
xlim([-5, 5]);
ylabel('Magnitude (V/m)');
title('FFT Plot of the Magnitude of Inputs and Outputs');
legend;
grid on;

% FFT Phase plot
figure('Name', 'FFT Phase');
plot(omega_THz, unwrap(angle(fftInputL)), 'c', 'DisplayName', 'Input', 'LineWidth', 2);
hold on;
plot(omega_THz, unwrap(angle(fftOutputR)), 'b', 'DisplayName', 'Right Output', 'LineWidth', 2);
plot(omega_THz, unwrap(angle(fftOutputL)), 'r', 'DisplayName', 'Left Output', 'LineWidth', 2);
hold off;
xlabel('Frequency (THz)');
ylabel('Phase (rad)');
title('FFT Phase of Input and Outputs');
legend;
grid on;
