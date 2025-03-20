% Name: Ayooluwa Owolabi and Student#: 101237575
% ELEC 4700 Winter 2025 semester
% Milestone 8:Lasers

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
c_eps_0 = 8.8542149e-12; % vacuum permittivity (F/m)
c_eps_0_cm = c_eps_0/100; % Vacuum permittivity measured in F/cm
c_mu_0 = 1/c_eps_0/c_c^2; % Permittivity of free space
c_q = 1.60217653e-19; % The value of an elementary charge
c_hb = 1.05457266913e-34;     % Dirac Constant
c_h = c_hb*2*pi; % Planck constant

% Input pulse parameters (Modulated Gaussian)
InputParasL.E0= 0*1e5;  % Amplitude
InputParasL.we= 0;    % Angular Frequency
InputParasL.t0= 30e-13;  % Time offset
InputParasL.wg= 10e-13;  % Gaussian width
InputParasL.phi= 0; % phase
InputParasL.rep = 0.5e-9;
InputParasR = 0;

% More constants definition
n_g = 3.5; % Group index of the medium
vg = c_c/n_g*1e2; % TWM Group velocity (cm/s)
L = 100e-4; %L = 1000e-6*1e2; % Medium length (cm)
Lambda = 1550e-9; % Wavelength of light
Nz = 51; % Spatial points
dz = L/(Nz-1); % Spatial step
dt = dz/vg; % Time step
Nt = floor(1000*Nz); % Time steps
f0 = c_c / Lambda;
Ntr = 1e18;
Nlim = [0, 3*Ntr];

% Define and Initialize the simulation parameters
plotN = 100; % How often the plots are updated (output speed)
XL = [0, L];
YL = [-4e6, 4e6]; % Range of the y-axis
fsync = dt*vg/dz; % Field sync factor
tmax = Nt*dt; % Max simulation time
t_L = dt*Nz; % time to travel the lenght of the medium

% Define spatial grid
z = linspace(0, L, Nz).';

% time = linspace(0, Nt*dt, Nt);
% Initialize input and output signals
time = nan(1, Nt);
InputL = nan(1, Nt);
InputR = nan(1, Nt);
OutputL = nan(1, Nt);
OutputR = nan(1, Nt);

% Forward and backward polarization of the electric field
Ef = zeros(size(z)); % forward E-field
Er = zeros(size(z)); % backward E-field
Pf = zeros(size(z)); % forward polarization
Pr = zeros(size(z)); % backward polarization

% Updated Polarization and propagation of the electric field
Efp = Ef;
Erp = Er;
Pfp = Pf;
Prp = Pr;

Ef1 = @SourceFct;
ErN = @SourceFct;

t = 0; % initial time to 0
time(1) = t;
InputL(1) = Ef1(t, InputParasL);
InputR(1) = ErN(t, InputParasR);
OutputR(1) = Ef(Nz);
OutputL(1) = Er(1);
Ef(1) = InputL(1);
Er(Nz) = InputR(1);

% Initialize variables
N = ones(size(z)) * Ntr;   % Initialize carrier density to Ntr across z
Nave(1) = mean(N);   % Track the average carrier density

% Define carrier injection parameters
gain = vg * 2.5e-16;   % Gain coefficient
eVol = 1.5e-10 * c_q;  % Effective volume of interaction
Ion = 0.05e-9;      % Injection current turn-on time
Ioff = 3000000e-9;  % Injection current turn-off time
I_off = 0.024;   % Carrier injection rate when source is off
I_on = 0.2;      % Carrier injection rate when source is on
taun = 1e-9;     % Carrier lifetime

% Compute the photon-to-carrier interaction scaling factor
Zg = sqrt(c_mu_0 / c_eps_0) / n_g;
EtoP = 1 / (Zg * f0 * vg * 1e-2 * c_hb);   % Convert field energy to photon density
alpha = 0;

% Stimulated emission couple
% beta_r = 0;
% gain_z = gain.*(N-Ntr)./vg;
% beta_i = (gain_z - alpha)./2;
% beta = beta_r + 1i*beta_i;
% exp_det = exp(-1i*dz*beta);
beta_spe = 0.3e-5;
gamma = 1.0;
SPE = 7;

% Reflection coefficients
RL = 0.5;
RR = 0.5;

% Define kappa (grating strength)
kappa0 = 0;
kappa = kappa0 * ones(size(z));
kappa(z < 0.3*L | z > 0.7*L) = 0; % Grating in central region

% Defining extra parameters for dispersion (gain/loss)
g_fwhm = 3.53e+012/10;
LGamma = g_fwhm*2*pi;
Lw0 = 0.0;
LGain = 0.015; %LGain = 0.01;

subplot(3,2,1)
plot(z*10000,real(Ef),'r'); % Real part of forward field
hold on
xlim(XL * 1e4);
ylim("auto");
xlabel('z(\mum)')
ylabel('E_f')
legend('\Re', '\Im');
hold off

subplot(3,2,2)
plot(z*10000,N, 'b');
hold on
xlim([0, L*10000]);
ylim(Nlim);
xlabel('z(\mum)')
ylabel('N')
hold off

subplot(3,2,[3,4])
plot(time*1e12, Nave, 'b');
xlim([0, Nt*dt*1e12])
ylim(Nlim)
xlabel('time(ps)')
ylabel('Nave')

subplot(3,2, [5,6])
plot(time*1e12,real(InputL),'r'); hold on
plot(time*1e12,real(OutputR), 'r--'); % Real part of input/output fields
plot(time*1e12,real(InputR),'b'); hold on
plot(time*1e12,real(OutputL),'b--');
xlim([0,Nt*dt*1e12])
ylim(YL)
xlabel('time(ps)')
ylabel('0')
hold off
legend('Left Input', 'Right Output', 'Right Input', 'Left Output', 'Location', 'east')
hold off

% Define loop for the time stepping of the field propagation
for i = 2:Nt
    t = dt * (i-1);
    time(i) = t;
    InputL(i) = Ef1(t, InputParasL);
    InputR(i) = ErN(t, 0);

    % Define Pf and Pr boundary conditions
    Pf(1) = 0;
    Pf(Nz) = 0;
    Pr(1) = 0;
    Pr(Nz) = 0;
    Cw0 = -LGamma + 1i*Lw0;

    S = (abs(Ef).^2 + abs(Er).^2).*EtoP*1e-6;
    if t < Ion || t > Ioff
        I_injv = I_off;
    else
        I_injv = I_on;
    end
    Stim = gain.*(N-Ntr).*S;
    N = (N + dt*(I_injv/eVol - Stim))./(1+ dt/taun);
    Nave(i) = mean(N(:));
    beta_r = 0;
    gain_z = gain.*(N-Ntr)./vg;
    beta_i = (gain_z - alpha)./2;
    beta = (beta_r + 1i*beta_i);
    exp_det = exp(-1i*dz*beta);

    A = sqrt(gamma*beta_spe*c_hb*f0*L*1e-2/taun)/(2*Nz);
    if SPE > 0
        eTf = (randn(Nz,1)+1i*randn(Nz,1))*A;
        eTr = (randn(Nz,1)+1i*randn(Nz,1))*A;
    else
        eTf = (ones(Nz,1))*A;
        eTr = (ones(Nz,1))*A;
    end

    EsF = eTf*abs(SPE).*sqrt(N.*1e6);
    EsR = eTr*abs(SPE).*sqrt(N.*1e6);

    % Computing Tf and Tr using the finite difference approximation
    Tf = LGamma*Ef(1:Nz-2) + Cw0*Pfp(2:Nz-1) + LGamma*Efp(1:Nz-2);
    Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5*dt*Tf)./(1-0.5*dt*Cw0);
    Tr = LGamma*Er(3:Nz) + Cw0*Prp(2:Nz-1) + LGamma*Erp(3:Nz);
    Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5*dt*Tr)./(1-0.5*dt*Cw0);

    % Modified Ef and Er to implement the reflections
    Ef(1) = InputL(i) + RL*Er(1);
    Er(Nz) = InputR(i) + RR * Ef(Nz);

    % Update forward and backward fields with grating terms
    Ef(2:Nz) = fsync*exp_det(1:Nz-1).*Ef(1:Nz-1) + 1i * dz * kappa(2:Nz) .* Er(2:Nz);
    Er(1:Nz-1) = fsync*exp_det(2:Nz).*Er(2:Nz) + 1i * dz * kappa(1:Nz-1) .* Ef(1:Nz-1);

    Ef(2:Nz-1) = Ef(2:Nz-1) - LGain * (Ef(2:Nz-1)- Pf(2:Nz-1));
    Er(2:Nz-1) = Er(2:Nz-1) - LGain * (Er(2:Nz-1)- Pr(2:Nz-1));

    % The modified Output implementing the mirrors
    OutputR(i) = Ef(Nz)*(1-RR);
    OutputL(i) = Er(1)*(1-RL);

    Ef = Ef + EsF;
    Er = Er + EsR;
    Efp = Ef;
    Erp = Er;
    Pfp = Pf;
    Prp = Pr;

    if mod(i,plotN) == 0
        
        subplot(3,2,1)
        plot(z*10000,real(Ef),'r'); % Real part of forward field
        hold on
        xlim(XL * 1e4);
        ylim("auto");
        xlabel('z(\mum)')
        ylabel('E_f')
        legend('\Re', '\Im');
        hold off

        subplot(3,2,2)
        plot(z*10000,N, 'b');
        hold on
        xlim([0, L*10000]);
        ylim(Nlim);
        xlabel('z(\mum)')
        ylabel('N')
        hold off

        subplot(3,2,[3,4])
        plot(time(1:i)*1e12, Nave(1:i), 'b')
        xlim([0, Nt*dt*1e12])
        ylim(Nlim)
        xlabel('time(ps)')
        ylabel('Nave')

        subplot(3,2, [5,6])
        plot(time*1e12,real(InputL),'r'); hold on
        plot(time*1e12,real(OutputR), 'r--'); % Real part of input/output fields
        plot(time*1e12,real(InputR),'b'); hold on
        plot(time*1e12,real(OutputL),'b--');
        xlim([0,Nt*dt*1e12])
        ylim(YL)
        xlabel('time(ps)')
        ylabel('0')
        hold off
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output', 'Location', 'east')
        hold off
        pause(0.01) % The pause considered to allow for the plotting animation
    end
end

% Compute and plot spectral content
fftOutputR = fftshift(fft(OutputR));
% fftOutputL = fftshift(fft(OutputL));
fftInputL = fftshift(fft(InputL));
omega = fftshift(wspace_M7(time));
figure('Name','Omega')
subplot(3,1,1)
plot(time*1e12,real(OutputR),'r'); hold on
plot(time*1e12,real(InputL),'g'); hold on  
legend('Output','Input')
xlabel('time(ps)')
ylabel('Output')
subplot(3,1,2)
plot(omega/2/pi*1E-12,20*log10(abs(fftOutputR))); hold on
plot(omega/2/pi*1E-12,20*log10(abs(fftInputL))); hold on
legend('Output','Input')
xlabel('GHz')
ylabel('20 log_{10} |E|')
xlim([-5 5]);
subplot(3,1,3)
plot(omega/2/pi*1E-12,unwrap(angle(fftOutputR))); hold on
plot(omega/2/pi*1E-12,unwrap(angle(fftInputL)))
xlabel('GHz')
ylabel('phase (E)')
xlim([-5 5]);
