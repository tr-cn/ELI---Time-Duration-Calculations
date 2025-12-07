clear; close all; clc;
%% Load spectrum: columns = [lambda_nm, I_lambda]
load('data.mat');   % <-- change filename as needed
lambda_nm = wavelengths_nm;               % wavelength in nm
I_lambda  = abs(y_pixels-max(y_pixels));               % spectral intensity (arb. units)

% Ensure column vectors and sort by wavelength
[lambda_nm, idx] = sort(lambda_nm(:));
I_lambda         = I_lambda(idx);

%% Convert wavelength axis to angular frequency axis
c = 299792458;
lambda_m = lambda_nm*1e-9;
omega   = 2*pi*c./lambda_m;          % non-uniform, descending

% Uniform omega grid
Nspec      = numel(omega);
omega_min  = min(omega);
omega_max  = max(omega);
omega_u    = linspace(omega_min, omega_max, Nspec).';
I_omega    = interp1(omega, I_lambda, omega_u, 'linear', 0);

%% Zero padding in frequency domain
ZP_factor  = 200;                      % interpolation factor in time
Nfft       = ZP_factor * Nspec;

% Centered spectrum, then pad with zeros
E_omega0   = sqrt(I_omega);         % flat phase
% E_omega0   = fftshift(E_omega0);

pad_total  = Nfft - Nspec;
pad_left   = floor(pad_total/2);
pad_right  = pad_total - pad_left;

E_omega_pad = [zeros(pad_left,1); E_omega0; zeros(pad_right,1)];

% New uniform omega grid matching padded array
domega     = omega_u(2) - omega_u(1);
omega_u_pad = ( (-(Nfft/2):(Nfft/2-1)).' ) * domega + mean(omega_u);

%% IFFT to time domain
% E_t = ifft(ifftshift(E_omega_pad)) * Nfft * domega / (2*pi);
E_t = ifft(E_omega_pad) * Nfft * domega / (2*pi);
% Time axis
T  = 2*pi/domega;                   % unchanged total window
dt = T/Nfft;                        % smaller time step due to padding
t  = (-Nfft/2:Nfft/2-1).' * dt;

%% Intensity and FWHM
I_t = abs(fftshift(E_t)).^2;

Imax = max(I_t);
half = Imax/2;
above = find(I_t >= half);
t_FWHM = t(above(end)) - t(above(1));

fprintf('Transform-limited pulse duration (FWHM) = %.3f fs\n', t_FWHM/1e-15);

%% Plots
figure;
subplot(2,1,1);
plot(lambda_nm, I_lambda);hold on;
% plot(I_omega)
xlabel('\lambda (nm)'); ylabel('I(\lambda)');
title('Input spectrum');

subplot(2,1,2);
plot(t*1e15, I_t);
xlabel('t (fs)'); ylabel('I(t)');
title(sprintf('Time-domain pulse, FWHM = %.2f fs', t_FWHM*1e15));