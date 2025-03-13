% Name: Ayooluwa Owolabi and Student#: 101237575
% ELEC 4700 Winter 2025 semester
% Milestone 6: Carrier Equation for Optical Amplifiers/ Lasers

clc; 
clear;
close all;

% Set default plot appearance and visualization settings
set(0, 'defaultaxesfontsize', 20) % Sets the default font size of the axes
set(0, 'DefaultFigureWindowStyle', 'docked') % Sets the default Figure Window style to be docked
set(0, 'DefaultLineLineWidth', 2) % Sets the Default Line Width; i.e. from 0 to 2
set(0, 'DefaultaxesLineWidth', 2) % Sets the Default Line Width of the axes

% Initializing Fundamental constants 
c_c = 299792458; % TWM speed of light (m/s)
c_eps_0 = 8.8542149e-12; % vacuum permittivity (F/m)
c_eps_0_cm = c_eps_0/100; % Vacuum permittivity measured in F/cm
c_mu_0 = 1/c_eps_0/c_c^2; % Permittivity of free space
c_q = 1.60217653e-19; % The value of an elementary charge
c_hb = 1.05457266913e-34;     % Dirac Constant
c_h = c_hb*2*pi; % Planck constant 
 
% Input pulse parameters (Modulated Gaussian)
InputParasL.E0 = 0*9e5;  % Amplitude (adjusted to match the scale in the image, e.g., ×10^6)
InputParasL.we = 200e-12;    % Angular Frequency
InputParasL.t0 = 10e-12;  % Time offset (in seconds, adjusted for ps scale)
InputParasL.wg = 1e-12;  % Pulse width (in seconds, adjusted for ps scale)
InputParasL.phi = 0; % phase
InputParasL.type = 'sinc'; % Initialize the new pulse type
InputParasL.rep = 500e-12; 
InputParasR = 0; % Right input

% Waveguide and optical parameters
n_g = 3.5; % Group index of the medium
vg = c_c/n_g*1e2; % TWM Group velocity (cm/s)
L = 1000e-6*1e2; % Medium length (cm)
Lambda = 1550e-9; % Wavelength of light
f0 = c_c / Lambda;
Nz = 100; % Spatial points
dz = L/(Nz-1); % Spatial step
dt = dz/vg; % Time step
Nt = floor(400*Nz); % Time steps

% Define spatial grid in cm, convert to µm for plotting
z = linspace(0, L, Nz).';

% Defining extra parameters (unchanged from original)
g_fwhm = 3.53e+012/10;
LGamma = g_fwhm*2*pi;
Lw0 = 0.0;
LGain = 0.05;

% Define kappa (grating strength)
kappa0 = 0;
kappa = kappa0 * ones(size(z));
kappa(z < 0.3 * L | z > 0.7 * L) = 0; % Grating in central region

% Plot kappa vs z
figure;
plot(z, kappa, 'k', 'LineWidth', 2);
xlabel('z (cm)');
ylabel('\kappa');
title('Grating Profile');

% Define carrier injection parameters
Ntr = 1e18;   % Carrier density threshold
gain = vg * 2.5e-16;   % Gain coefficient
eVol = 1.5e-10 * c_q;  % Effective volume of interaction
Ion = 0.25e-9;   % Injection current turn-on time
Ioff = 3e-9;     % Injection current turn-off time
I_off = 0.024;   % Carrier injection rate when source is off
I_on = 0.1;      % Carrier injection rate when source is on
taun = 1e-9;     % Carrier lifetime

% Compute the photon-to-carrier interaction scaling factor
Zg = sqrt(c_mu_0 / c_eps_0) / n_g;
EtoP = 1 / (Zg * f0 * vg * 1e-2 * c_hb);   % Convert field energy to photon density
alpha = 0;

% ===============================================================
% Initialize variables 
N = ones(size(z)) * Ntr;   % Initialize carrier density to Ntr across z
Nave = nan(1, Nt);
Nave(1) = mean(N);   % Track the average carrier density

% Initialize input and output signals
time = nan(1, Nt);
InputL = zeros(1, Nt);
InputR = zeros(1, Nt);
OutputL = zeros(1, Nt);
OutputR = zeros(1, Nt);

Ef = ones(size(z)) * 1e3; % Forward propagating electric field
Er = ones(size(z)) * 1e5; % Reverse propagating electric field
fsync = dt*vg/dz; % Field sync factor

% Photon density S(z)
S = zeros(size(z)); % Initialize photon density

figure(1);
time = linspace(0, Nt * dt, Nt);

% Define loop for the time stepping of the field propagation
for i = 2:Nt
    time_vec(i) = dt * (i-1);  % Update time at each step (i is now explicitly an integer)
    t = time_vec(i); 
    InputL(i) = SourceFct_M6(t, InputParasL);
    InputR(i) = 0;

    % The carrier equation (Photon density)
    S = (abs(Ef).^2 + abs(Er).^2).*EtoP*1e-6;

    if (t < Ion || t > Ioff)
        I_injv = I_off;
    else
        I_injv = I_on;
    end

    % Modified and updated carrier density 
    Stim = gain.*(N - Ntr).*S;
    N = (N + dt*(I_injv/ eVol - Stim))./(1 + dt / taun);
    N(N < 0) = 0;
    Nave(i) = mean(N);


    % Forward and backward propagation of the fields 
    Ef(2:Nz) = fsync * Ef(1:Nz-1) + InputL(i);
    Er(1:Nz-1) = fsync * Er(2:Nz) + InputR(i);

    % The modified Output 
    OutputR(i) = Ef(Nz);
    OutputL(i) = Er(1);

    % Periodic plotting of the field
    if mod(i, 1000) == 0

        % Define common properties for plots
        lineWidth = 2;
        gridOn = 'on';
        xLabelFontSize = 12;
        yLabelFontSize = 12;
        titleFontSize = 14;

        % Plot 1: Carrier Density Profile
        subplot(2, 2, 1);
        plot(z * 1e4, N, 'r', 'LineWidth', lineWidth);
        xlabel('z (\mum)', 'FontSize', xLabelFontSize);
        ylabel('N(z)', 'FontSize', yLabelFontSize);
        title('Carrier Density N(z) Over Time', 'FontSize', titleFontSize);
        grid(gridOn);
        xlim([0, max(z) * 1e4]);
        ylim([0, 5e18]);  % Adjusted limits to show increase and decrease

        % Plot 2: Average Carrier Density Over Time
        subplot(2, 2, 2);
        plot(time(1:i) * 1e12, Nave(1:i), 'b', 'LineWidth', lineWidth);
        xlabel('Time (ps)', 'FontSize', xLabelFontSize);
        ylabel('N_{ave}', 'FontSize', yLabelFontSize);
        title('Average Carrier Density Over Time', 'FontSize', titleFontSize);
        grid(gridOn);
        xlim([0, 5000]);
        ylim([0, 5e18]);  % Adjusted limits to show rise and fall

        % Plot 3: Photon Density S(z)
        subplot(2, 2, 3);
        plot(z * 1e4, S, 'm', 'LineWidth', lineWidth);
        xlabel('z (\mum)', 'FontSize', xLabelFontSize);
        ylabel('S(z)', 'FontSize', yLabelFontSize);
        title('Photon Density S(z) Over Time', 'FontSize', titleFontSize);
        grid(gridOn);
        xlim([0, max(z) * 1e4]);
        ylim([0, max(S) + 1e6]);

        % Plot 4: Overlay of Input/Output Signals
        subplot(2, 2, 4);
        plot(time(1:i) * 1e12, real(InputL(1:i)), 'r', 'LineWidth', lineWidth); 
        hold on;
        plot(time(1:i) * 1e12, real(OutputR(1:i)), 'g', 'LineWidth', lineWidth);
        plot(time(1:i) * 1e12, real(InputR(1:i)), 'b', 'LineWidth', lineWidth);
        plot(time(1:i) * 1e12, real(OutputL(1:i)), 'm', 'LineWidth', lineWidth);
        xlabel('Time (ps)', 'FontSize', xLabelFontSize);
        ylabel('O', 'FontSize', yLabelFontSize);
        title('Input and Output Signals', 'FontSize', titleFontSize);
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output', 'Location', 'best');
        grid(gridOn);
        xlim([0, 5000]);
        ylim([-1e6, 1e7]);
        hold off

        pause(0.01); % Allow for animation
    end
end

% =========================================================================
% % Define and initialize simulation parameters
% plotN = 10; % How often the plots are updated
% XL = [0, L];
% YL = [-1 * InputParasL.E0, 2 * InputParasL.E0]; % Range of the y-axis
% fsync = dt * vg / dz; % Field sync factor
% tmax = Nt * dt; % Max simulation time
% t_L = dt * Nz; % Time to travel the length of the medium
% % Define Gaussian-modulated input source
% SourceFct = @(t, InputParas) InputParas.E0 * exp(-((t-InputParas.t0).^2) / (2 * InputParas.wg^2)) .* exp(1j * (InputParas.we * t + InputParas.phi));
% % Reflection coefficients
% RL = 0i;
% RR = 0i;
% % Adding more input parameters and simulation parameters (Gain and Detuning
% % factors)
% beta_r = 0; % Gain Factor
% beta_i = 0; % Detuning factor
% % Modify the plotting section to match the image
% plotN = 10; % How often the plots are updated
% XL = [0, L*1e4]; % z in µm
% YL = [-1*InputParasL.E0, 2*InputParasL.E0]; % Range of y-axis for electric fields
% tmax = Nt * dt; % Max simulation time
% t_L = dt * Nz; % time to travel the length of the medium
% % Initialize the dispersion terms
% Pf = zeros(size(z));
% Pr = zeros(size(z));
% % The previous field terms
% Efp = Ef;
% Erp = Er;
% Pfp = Pf;
% Prp = Pr;
% % Modified Ef and Modified Er
% beta = ones(size(z))*(beta_r+1i*beta_i);
% exp_det = exp(-1i*dz*beta);
% =========================================================================

% =========================================================================
% From the for loop
    % Stores the previous Ef values
    % % Modified Ef and Er to implement the reflections (unchanged)
    % Ef(1) = InputL(i) + RL*Er(1);
    % Er(Nz) = RR * Ef(Nz);
    % Ef_prev = Ef;
    % Er_prev = Er;
    % % Initialize polarization variables 
    % P = zeros(Nz, 1);
    % omega0 = 2 * pi * c_c / Lambda; % Resonant frequency
    % gamma = 1e12; % Damping factor
    % chi = 1e-3; % Susceptibility
    % % Update polarization using Backward Euler method 
    % P = (P + dt * chi .* Ef) ./ (1 + dt * (gamma - 1i * omega0));    
    % 
    % % Define Pf and Pr boundary conditions 
    % Pf(1) = 0;
    % Pf(Nz) = 0;
    % Pr(1) = 0;
    % Pr(Nz) = 0;
    % Cw0 = -LGamma + 1i*Lw0;
    % 
    % % Computing Tf and Tr using the finite difference approximation 
    % Tf = LGamma*Ef(1:Nz-2) + Cw0*Pfp(2:Nz-1) + LGamma*Efp(1:Nz-2);
    % Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5*dt*Tf)./(1-0.5*dt*Cw0);
    % Tr = LGamma*Er(3:Nz) + Cw0*Prp(2:Nz-1) + LGamma*Erp(3:Nz);
    % Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5*dt*Tr)./(1-0.5*dt*Cw0);
    % 
    % Ef(2:Nz-1) = Ef(2:Nz-1) - LGain * (Ef(2:Nz-1)- Pf(2:Nz-1));
    % Er(2:Nz-1) = Er(2:Nz-1) - LGain * (Er(2:Nz-1)- Pr(2:Nz-1)); 
% ========================================================================= 
    