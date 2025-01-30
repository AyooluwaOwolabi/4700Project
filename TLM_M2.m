% Name: Ayooluwa OWOLABI and Std#: 101237575 
% ELEC 4700A  Project Milestone1
%Set default plot appearance and visualization settings
set(0, 'defaultaxesfontsize',20) % Sets the default font size of the axes
set(0, 'DefaultFigureWindowStyle','docked') % Sets the default Figure Window style to be docked 
set(0, 'DefaultLineLineWidth',2) % Sets the Default Line Width; i.e. from 0 to 2 
set(0, 'DefaultaxesLineWidth',2) % Sets the Default Line Width of the axes

% set(0, 'DefaultFigureWindowStyle','docked') - Duplicated code 

% Initializing the fundamental constants 
c_c = 299792458; % TWM speed of light (m/s)
c_eps_0 = 8.8542149e-12; % vacuum permittivity (F/m)
c_eps_0_cm = c_eps_0/100; % Vacuum permittivity measured in F/cm
c_mu_0 = 1/c_eps_0/c_c^2; % Permittivity of free space
c_q = 1.60217653e-19; % The value of an elementary charge
c_hb = 1.05457266913e-34;     % Dirac Constant
c_h = c_hb*2*pi; % Planck constant 

% % Define and Initialize the input parameters for the left and right sources
% InputParasL.E0= 1e5;  % Amplitude
% InputParasL.we= 0;    % Angular Frequency
% InputParasL.t0= 2e-12;  % Time offset 
% InputParasL.wg= 5e-13;  % Gaussian width 
% InputParasL.phi= 0; % phase 
% InputParasR= 0; % No input is implemented for the right source  

% Define and Initialize the input parameters for the left and right sources
% (Modulated Gaussian)
InputParasL.E0= 1e5;  % Amplitude
InputParasL.we= 5e12;    % Angular Frequency
InputParasL.t0= 2e-12;  % Time offset 
InputParasL.wg= 5e-13;  % Gaussian width 
InputParasL.phi= 0; % phase 
InputParasR= 0; % No input is implemented for the right source 

% The source function with modulation
function E = SourceFct_M2(t, InputParas)
E = InputParas.E0 * exp(-((t-InputParas.t0.^2) / (2 * inputparas.wg^2)).* exp(1j * (Inputparas.we * t + InputParas.phi)));
end 

% Adding the reflection coefficients for mirrors 
RL = 0.9i; % The left mirror reflection coefficient 
RR = 0.9i; % The right mirror reflection coefficient 

% Adding more input parameters for the gain/loss
beta_r = 80;
beta_i =8;

% Define material properties and simulation parameters
n_g = 3.5; % Group index of the medium
vg = c_c/n_g*1e2; % TWM group velocity (cm/s)
Lambda = 1550e-9; % Wavelength of light (m)

% Define and Initialize the simulation parameters
% plotN = 10; % How frequent the plots are updated; Output speed 
plotN = 10; % How often the plots are updated has been increased to 50; making the resulting plot output to appear faster
L = 1000e-6*1e2; % Length of the medium
XL = [0,L]; % Range for the z-axis
YL = [-2*InputParasL.E0, 5*InputParasL.E0]; % Range for the y-axis

Nz = 500; % Number of points of the z-axis
dz = L/(Nz-1); % The spatial step size
dt = dz/vg; % Time step size 
fsync = dt*vg/dz; % Field sync factor 

Nt = floor(2*Nz); % Number of Time steps 
tmax = Nt*dt; % Max simulation time
t_L = dt*Nz;   % time to travel the length of the medium

% Initialize the spatial (z-axis) and temporal arrays 
z = linspace(0,L,Nz).';  %Nz points, Nz-1 segments 
time = nan(1,Nt); 
InputL = nan(1,Nt);
InputR = nan(1,Nt);
OutputL = nan(1,Nt);
OutputR = nan(1,Nt);

Ef = zeros(size(z)); % Forward propagating electric field 
Er = zeros(size(z)); % Reverse propagating electric field 

% Modified Ef and Modified Er 
beta = ones(size(z))*(beta_r+1i*beta_i);
exp_det = exp(-1i*dz*beta);

Ef1 = @SourceFct; % Left source function
ErN = @SourceFct; % Right source function

% The intial time step
t = 0; 
time(1) = t;

% Set initial conditions for the fields 
InputL(1) = Ef1(t,InputParasL); % Left input 
InputR(1) = ErN(t,InputParasR); % Right input 
OutputR(1) = Ef(Nz); % Right output 
OutputL(1) = Er(1); % Left output 
Ef(1) = InputL(1); % Forward field 
Er(Nz) = InputR(1); % Reverse field

% Plot initial fields 
figure('name','Fields')
subplot(3,1,1)
plot(z*10000,real(Ef),'r'); % Real part of forward field 
hold off
xlabel('z(\mum)')
ylabel('E_f')
subplot(3,1,2)
plot(z*10000,real(Er), 'b'); % Real part of the reverse field 
xlabel('z(\mum)')
ylabel('E_r')
hold off
subplot(3,1,3)
plot(time*1e12,real(InputL),'r'); hold on
plot(time*1e12,real(OutputR), 'r--'); % Real part of input/output fields
plot(time*1e12,real(InputR),'b'); hold on
plot(time*1e12,real(OutputL),'b--');
xlabel('time(ps)')
ylabel('E')

hold off

% Define loop for the time stepping of the field propagation 
for i = 2:Nt
    % Modified Ef and Er to implement the reflections 
    Ef(1) = InputL(i) + RL*Er(1);
    Er(Nz) = InputR(i) + RR*Ef(Nz);
    
    % Modified Output for both the right and left of the mirror 
    OutputR(i) = Ef(Nz)*(1-RR);
    OutputL(i) = Er(1)*(1-RL);

    % The time is updated 
    t = dt*(i-1);
    time(i) = t;

    InputL(i) = Ef1(t,InputParasL); % Left input field
    InputR(i) = ErN(t, 0); % Right input field

    % The boundary conditions
    Ef(1) = InputL(i);
    Er(Nz) = InputR(i);

    % Modifying Ef and Er to implement mirrors 
    Ef(1) = InputL(i) + RL*Er(1);
    Er(Nz) = InputR(i) + RR*Ef(Nz);

    % The forward and backward propagations of the field
    % Ef(2:Nz) = fsync*Ef(1:Nz-1); % The forward propagation for the field  
    % Er(1:Nz-1) = fsync*Er(2:Nz); % The reverse propagation (backwards) for the field 

    % Updating the forward and backward propagations of the field 
    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz-1);
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz);

    % Updating the output so it reflect on the field 
    OutputR(i) = Ef(Nz);
    OutputL(i) = Er(1);

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


% Fourier Transform of OutputR 
fftOutput = fftshift(fft(OutputR));
omega = fftshift(wspace_M2(time));

% Plot the time-domain signal
figure; 
plot(time * 1e12, real(OutputR), 'r');
xlabel('Time (ps)');
ylabel('Output Signal');
title('Time Domain Output vs Time');

figure;
tiledlayout(2,1);

% Plot magnitude of FFT (Top plot)
nexttile;
plot(omega, abs(fftOutput), 'm');
xlabel('Frequency (rad/s)');
ylabel('|FFT(OutputR)|');
title('Magnitude of FFT Output');

% Plot phase of FFT (using unwrap) (Bottom plot)
nexttile;
plot(omega, unwrap(angle(fftOutput)), 'r');
xlabel('Frequency (rad/s)');
ylabel('Phase (radians)');
title('Phase of FFT Output');