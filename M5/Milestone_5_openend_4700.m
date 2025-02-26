% Name: Ayooluwa Owolabi and Student#: 101237575
% Milestone 5: Open-ended Milestone for First Project Report 

clc; 
clear;
close all;

% Set default plot appearance and visualization settings
set(0, 'defaultaxesfontsize',20) % Sets the default font size of the axes
set(0, 'DefaultFigureWindowStyle','docked') % Sets the default Figure Window style to be docked
set(0, 'DefaultLineLineWidth',2) % Sets the Default Line Width; i.e. from 0 to 2
set(0, 'DefaultaxesLineWidth',2) % Sets the Default Line Width of the axes

% Initializing Fundamental constants
c_c = 299792458; % Speed of light (m/s)
n_g = 3.5; % Group index of the medium
vg = c_c/n_g*1e2; % TWM Group velocity (cm/s)
L = 1000e-6*1e2; % Medium length (cm)
Lambda = 1550e-9; % Wavelength of light
Nz = 500; % Spatial points
dz = L/(Nz-1); % Spatial step
dt = dz/vg; % Time step
Nt = floor(2*Nz); % Time steps
time = linspace(0, Nt*dt, Nt);

% Defining extra parameters
g_fwhm = 3.53e+012/10;
LGamma = g_fwhm*2*pi;
Lw0 = 0.0;
LGain = 0.01;

% Define a range of time step values to vary fsync
dt_val = [dz/(vg * 0.5), dz/vg, dz/(vg * 1.5)]; % Differing values from the accurate fsync value
% dt_val = linspace(dz/(vg*0.5), dz/(vg*2), 5); % More variations of dt % More variations of dt 


% Input pulse parameters (Modulated Gaussian)
InputParasL.E0= 1e5;  % Amplitude
InputParasL.we= 0;    % Angular Frequency
InputParasL.t0= 2e-12;  % Time offset
InputParasL.wg= 5e-13;  % Gaussian width
InputParasL.phi= 0; % phase
InputParasR = 0;

% Define spatial grid
z = linspace(0, L, Nz).';

% ============== The old code for implementation of grating ===========
% % Define kappa (grating strength)
% kappa0 = 0;
% % kappa = kappa0 * ones(size(z));
% kappa(z < 0.3*L | z > 0.7*L) = 0; % Grating in central region

% % Plot kappa vs z
% figure;
% plot(z, kappa, 'k', 'LineWidth', 2);
% xlabel('z (cm)');z
% ylabel('\kappa');
% title('Grating Profile');
% ===============================================================

% ====== Define the structured grating for bit encoding =========

% Define structured grating for bit encoding
bit_sequences = {[1 0 1 1 0], [1 1 0 1 0 1], [0 1 1 0 1]}; % Multiple bit patterns
for seq_idx = 1:length(bit_sequences)
    bit_pattern = bit_sequences{seq_idx};
    bit_length = length(bit_pattern);
    kappa = zeros(Nz,1); % The new kappa used for the bit encoding
    bit_width = floor(Nz / bit_length);

    for i = 1:bit_length
        if bit_pattern(i) == 1
            kappa((i-1)*bit_width + (1:bit_width)) = 0.05; % Grating strength for '1'
        end
    end
end 
% ===============================================================

% Initialize fields
Ef = zeros(size(z));
Er = zeros(size(z));

% Define Gaussian-modulated input source
SourceFct = @(t, InputParas) InputParas.E0 * exp(-((t-InputParas.t0).^2) / (2 * InputParas.wg^2)) .* exp(1j * (InputParas.we * t + InputParas.phi));

% Reflection coefficients
RL = 0i;
RR = 0i;

% Adding more input parameters and simulation parameters (Gain and Detuning factors)
beta_r = 0; % Gain Factor
beta_i = 0; % Detuning factor

% Define and Initialize the simulation parameters
plotN = 10; % How often the plots are updated (output speed)
XL = [0, L];
YL = [-1*InputParasL.E0, 2*InputParasL.E0]; % Range of the y-axis
fsync = dt*vg/dz; % Field sync factor
tmax = Nt * dt; % Max simulation time
t_L = dt * Nz; % time to travel the length of the medium

% Initialize input and output signals
InputL = zeros(1, Nt);
InputR = zeros(1, Nt);
OutputL = zeros(1, Nt);
OutputR = zeros(1, Nt);

% Ef = zeros(size(z)); % Forward propagating electric field
% Er = zeros(size(z)); % Reverse propagating electric field

% Initialize the dispersion terms
Pf = zeros(size(z));
Pr = zeros(size(z));

% The previous field terms
Efp = Ef;
Erp = Er;
Pfp = Pf;
Prp = Pr;

% Modified Ef and Modified Er
beta = ones(size(z))*(beta_r+1i*beta_i);
exp_det = exp(-1i*dz*beta);

% Initialize time vector
time_vec = zeros(1, Nt);  % Allocate memory to avoid dynamic resizing

% Define loop for the time stepping of the field propagation
for i = 2:Nt
    time_vec(i) = dt * (i-1);  % Update time at each step
    t = time_vec(i);
    InputL(i) = SourceFct(t, InputParasL);

    % Modified Ef and Er to implement the reflections
    Ef(1) = InputL(i) + RL*Er(1);
    Er(Nz) = RR * Ef(Nz);

    % Stores the previous Ef values
    Ef_prev = Ef;
    Er_prev = Er;
    
    % Initialize polarization variables
    P = zeros(Nz, 1);
    omega0 = 2 * pi * c_c / Lambda; % Resonant frequency
    gamma = 1e12; % Damping factor
    chi = 1e-3; % Susceptibility

    % Update polarization using Backward Euler method (implicit method)
    P = (P + dt * chi .* Ef) ./ (1 + dt * (gamma - 1i * omega0));
    
    % ============ Before including the updated polarization scheme ==============
    % % Updated forward and backward fields with the grating terms implemented using bit encoding 
    % Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz-1) + 1i * dz * kappa(2:Nz) .* Er(2:Nz);
    % Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i * dz * kappa(1:Nz-1) .* Ef_prev(1:Nz-1);

    % Forward and backward propagation of the fields with grating terms implemented using bit encoding with improved boundary handling
    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz-1) + 1i * dz * kappa(2:Nz) .* Er(2:Nz) + P(2:Nz);
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i * dz * kappa(1:Nz-1) .* Ef_prev(1:Nz-1) + P(1:Nz-1);

    % Define Pf and Pr boundary conditions
    Pf(1) = 0;
    Pf(Nz) = 0;
    Pr(1) = 0;
    Pr(Nz) = 0;
    Cw0 = -LGamma + 1i*Lw0;

    % Computing Tf and Tr using the finite difference approximation
    Tf = LGamma*Ef(1:Nz-2) + Cw0*Pfp(2:Nz-1) + LGamma*Efp(1:Nz-2);
    Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5*dt*Tf)./(1-0.5*dt*Cw0);
    Tr = LGamma*Er(3:Nz) + Cw0*Prp(2:Nz-1) + LGamma*Erp(3:Nz);
    Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5*dt*Tr)./(1-0.5*dt*Cw0);

    Ef(2:Nz-1) = Ef(2:Nz-1) - LGain * (Ef(2:Nz-1)- Pf(2:Nz-1));
    Er(2:Nz-1) = Er(2:Nz-1) - LGain * (Er(2:Nz-1)- Pr(2:Nz-1));

    % The modified Output implementing the mirrors
    OutputR(i) = Ef(Nz)*(1-RR);
    OutputL(i) = Er(1)*(1-RL);
    
    % The periodic plotting of the field
    if mod(i,plotN) == 0
        subplot(3,1,1)
        plot(z*10000,real(Ef), 'r'); hold on
        plot(z*10000,imag(Ef), 'r--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_f')
        legend('\Re','\Im')
        hold off
        subplot(3,1,2)
        plot(z*10000,real(Er), 'b'); hold on
        plot(z*10000,imag(Er), 'b--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_r')
        legend('\Re', '\Im')

        hold off
        subplot(3,1,3);
        plot(time*1e12,real(InputL), 'r'); hold on
        plot(time*1e12,real(OutputR), 'g');
        plot(time*1e12,real(InputR), 'b');
        plot(time*1e12,real(OutputL), 'm');
        xlim([0,Nt*dt*1e12])
        ylim(YL)
        xlabel('time(ps)')
        ylabel('0')
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output', 'Location', 'east')
        hold off
        pause(0.01) % The pause considered to allow for the plotting animation
    end

    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;
end

% ============================ Grating Plots ==============================
% Plot Grating Configuration
figure;
subplot(3,1,1);
plot(linspace(0, L, Nz) * 1e4, real(kappa), 'm');
xlabel('z (\mum)');
ylabel('Re(\kappa)');
% title('Grating Configuration');
title(['Grating Configuration for Bit Sequence ', num2str(seq_idx)]);

% Pulse Propagation Plot
subplot(3,1,2);
plot(linspace(0, L, Nz) * 1e4, abs(Ef), 'm');
hold on;
plot(linspace(0, L, Nz) * 1e4, abs(Er), 'r');
xlabel('z (\mum)');
ylabel('|E_f|, |E_r|');
title('Pulse Propagation');

% % Overlay Pulse Generation for different dt values
% subplot(3,1,3);
% plot(time * 1e12, real(OutputR), 'c');
% hold on;
% legend({'\deltat = 30fs, \deltas = 75fs', '\deltat = 30fs, \deltas = 200fs', '\deltat = 75fs, \deltas = 300fs'});
% xlabel('time (ps)');
% ylabel('|e_{out}|');
% title('Pulse generation for different synchronization conditions');
% hold off;

% Overlay Pulse Generation for different dt values
subplot(3,1,3);
plot(time * 1e12, real(OutputR), 'c');
% legend({'\deltat = 30fs, \deltas = 75fs', '\deltat = 30fs, \deltas = 200fs', '\deltat = 75fs, \deltas = 300fs'});
legend(arrayfun(@(x) sprintf('\Deltat = %.1fft', x*1e15), dt_val, 'UniformOutput', false));
xlabel('Time (ps)');
ylabel('|E_{out}|');
title('Pulse Generation for Different Synchronization Conditions');
hold off;

% =========================================================================

% Compute and plot spectral content
fftOutputR = fftshift(fft(OutputR));
fftOutputL = fftshift(fft(OutputL));
fftInputL = fftshift(fft(InputL));
omega = fftshift(wspace_M4(time));

% Plot time-domain output
figure;
plot(time_vec * 1e12, real(OutputR), 'r');
xlabel('Time (ps)');
ylabel('Output Signal');
title('Time Domain Output vs Time');

% -------------------------------------------------------------------------
% This is to show the plot for the modified synchronized fsync
figure; 
hold on;
for idx = 1:length(dt_val)
    dt = dt_val(idx);
    time = linspace(0, Nt * dt, Nt);
    fsync = dt * vg / dz; % Modified synchronization factor
    beta = 2 * pi / Lambda; % Propagation constant 

    % Initialize polarization variables
    P = zeros(Nz, 1);
    omega0 = 2 * pi * c_c / Lambda; % Resonant frequency
    gamma = 1e12; % Damping factor
    chi = 1e-3; % Susceptibility

    % subplot(2,1,1) % 2 rows 1 colummn and first position
    subplot(length(dt_val), 2, 2*idx-1);
    plot(time_vec * 1e12, real(OutputR), 'b');
    xlabel('Time (ps)');
    ylabel('Output Signal');
    title(['Synchronization Condition: fsync = ' num2str(fsync)]);

    % Plotting the frequency content
    % subplot(2,1,2) % 2 rows 1 colummn and second position
    subplot(length(dt_val), 2, 2*idx);
    plot(omega, abs(fftOutputR), 'm');
    xlabel('Frequency (rad/s)');
    ylabel('Magnitude');
    title(['Frequency Content: fsync = ' num2str(fsync)]);
end

% Alternative: Simple Amplitude Modulation for Comparison
figure;
am_mod = SourceFct(time, InputParasL) .* cos(2 * pi * (c_c / Lambda) * time);
plot(time * 1e12, am_mod, 'k');
hold on;
plot(time * 1e12, real(OutputR), 'c');
legend({'Amplitude Modulation', 'Structured Grating Output'});
xlabel('Time (ps)');
ylabel('Signal Amplitude');
title('Comparison: Simple Modulation vs. Structured Grating');

% -------------------------------------------------------------------------

% =========== Plotting the required graphs on the same plot ================
figure;
plot(omega, abs(fftOutputR), 'r', 'LineWidth', 1.5); hold on;
plot(omega, abs(fftOutputL), 'b', 'LineWidth', 1.5);
plot(omega, abs(fftInputL), 'm', 'LineWidth', 1.5);
hold off;

xlabel('Frequency (rad/s)');
ylabel('Magnitude');
title('FFT Magnitude Comparison');
legend('FFT Output (Right) - High Pass', 'FFT Output (Left) - Low Pass', 'FFT Input (Gaussian)');
grid on;

% Display that the simulation is complete 
disp('Simulation complete. Check for pulse distortion due to the fsync variation.');
