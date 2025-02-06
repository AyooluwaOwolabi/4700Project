% ELEC 4700A - Project Milestone 3
% Name: Ayooluwa Owolabi and Student#: 101237575
% Adding a grating to the waveguide and analyzing wave propagation

% Set default plot appearance and visualization settings 
set(0, 'defaultaxesfontsize',20) % Sets the default font size of the axes
set(0, 'DefaultFigureWindowStyle','docked') % Sets the default Figure Window style to be docked 
set(0, 'DefaultLineLineWidth',2) % Sets the Default Line Width; i.e. from 0 to 2 
set(0, 'DefaultaxesLineWidth',2) % Sets the Default Line Width of the axes
% set(0, 'DefaultFigureWindowStyle','docked') - Duplicated code 


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

% Input pulse parameters (Modulated Gaussian)
InputParasL.E0= 1e5;  % Amplitude
InputParasL.we= 0;    % Angular Frequency
InputParasL.t0= 2e-12;  % Time offset 
InputParasL.wg= 5e-13;  % Gaussian width 
InputParasL.phi= 0; % phase
InputParasR = 0;

% Define spatial grid
z = linspace(0, L, Nz).';

% Define kappa (grating strength)
kappa0 = 100;
kappa = kappa0 * ones(size(z));
kappa(z < 0.3*L | z > 0.7*L) = 0; % Grating in central region

% Plot kappa vs z
figure;
plot(z, kappa, 'k', 'LineWidth', 2);
xlabel('z (cm)');
ylabel('\kappa');
title('Grating Profile');

% Initialize fields
Ef = zeros(size(z));
Er = zeros(size(z));

% Define Gaussian-modulated input source
SourceFct = @(t, InputParas) InputParas.E0 * exp(-((t-InputParas.t0).^2) / (2 * InputParas.wg^2)) .* exp(1j * (InputParas.we * t + InputParas.phi));

% Reflection coefficients
RL = 0i;
RR = 0i;

% Adding more input parameters and simulation parameters (Gain and Detuning factors)
beta_r = 0; % Gain Factor (γ)
beta_i = 0; % Detuning factor ​
% Define and Initialize the simulation parameters 
plotN = 50; % How often the plots are updated (output speed)
XL = [0, L];
YL = [-2*InputParasL.E0, 5*InputParasL.E0]; % Range of the y-axis 
fsync = dt*vg/dz; % Field sync factor
tmax = Nt*dt; % Max simulation time 
t_L = dt*Nz; % time to travel the lenght of the medium 

% Initialize input and output signals
InputL = zeros(1, Nt); 
InputR = zeros(1, Nt);
OutputL = zeros(1, Nt); 
OutputR = zeros(1, Nt);

% Ef = zeros(size(z)); % Forward propagating electric field 
% Er = zeros(size(z)); % Reverse propagating electric field 

% Modified Ef and Modified Er 
beta = ones(size(z))*(beta_r+1i*beta_i);
exp_det = exp(-1i*dz*beta);

% Ef1 = @SourceFct_M3; % Left source function
% ErN = @SourceFct_M3; % Right source function
% 
% % The intial time step
% t = 0; 
% time(1) = t;
% 
% % Set initial conditions for the fields 
% InputL(1) = Ef1(t,InputParasL); % Left input 
% InputR(1) = ErN(t,InputParasR); % Right input 
% OutputR(1) = Ef(Nz); % Right output 
% OutputL(1) = Er(1); % Left output 
% Ef(1) = InputL(1); % Forward field 
% Er(Nz) = InputR(1); % Reverse field
% 
% % Plot initial fields 
% figure('name','Fields')
% subplot(3,1,1)
% plot(z*10000,real(Ef),'r'); % Real part of forward field 
% hold off
% xlabel('z(\mum)')
% ylabel('E_f')
% subplot(3,1,2)
% plot(z*10000,real(Er), 'b'); % Real part of the reverse field 
% xlabel('z(\mum)')
% ylabel('E_r')
% hold off
% subplot(3,1,3)
% plot(time*1e12,real(InputL),'r'); hold on
% plot(time*1e12,real(OutputR), 'r--'); % Real part of input/output fields
% plot(time*1e12,real(InputR),'b'); hold on
% plot(time*1e12,real(OutputL),'b--');
% xlabel('time(ps)')
% ylabel('E')
% 
% hold off

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

    % % The time is updated 
    % t = dt*(i-1);
    % time(i) = t;

    % Stores the previous Ef values 
    Ef_prev = Ef;
    Er_prev = Er;
    
    % Update forward and backward fields with grating terms
    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz-1) + 1i * dz * kappa(2:Nz) .* Er(2:Nz);
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i * dz * kappa(1:Nz-1) .* Ef_prev(1:Nz-1);

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
end 

% Compute and plot spectral content
fftOutputR = fftshift(fft(OutputR));
fftOutputL = fftshift(fft(OutputL));
fftInputL = fftshift(fft(InputL));
omega = fftshift(wspace_M3(time));

% Plot time-domain output
figure;
plot(time_vec * 1e12, real(OutputR), 'r');
xlabel('Time (ps)'); 
ylabel('Output Signal');
title('Time Domain Output vs Time');

figure;
tiledlayout(3,1);

% Plot FFT Magnitude (Right Output)
% This gives a high pass filter (transmission) because the low frequencies are attenuated
% allowing the high frequencies to match and pass through
nexttile;
plot(omega, abs(fftOutputR), 'r');
xlabel('Frequency (rad/s)');
ylabel('|FFT(OutputR)|');
title('Magnitude of FFT Output (Right)');

% Plot FFT Magnitude (Left Output)
% This gives a low pass filter (reflection) because the high frequencies
% are attenuated allowing the low frequencies to match and pass through
% attenuated 
nexttile;
plot(omega, abs(fftOutputL), 'b');
xlabel('Frequency (rad/s)');
ylabel('|FFT(OutputL)|');
title('Magnitude of FFT Output (Left)');

% Plot FFt Magnitude (Gaussian)
nexttile; 
plot(omega, abs(fftInputL), 'm');
xlabel('Frequency (rad/s)');
ylabel('|InputL|');
title('Magnitude of InputL')
