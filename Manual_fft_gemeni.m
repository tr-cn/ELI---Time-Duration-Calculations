clear; close all; clc;
% --- 1. Load and Convert Data ---
load("Real_data.mat")

c = 299792458; % Speed of light (m/s)

% Your Input Data (Wavelength in nm, Spectrum is REAL Intensity)
wavelength_nm = wavelengths_nm; % Example: Your input vector in nm
I_spectrum = abs(y_pixels - max(y_pixels));

% Convert to Angular Frequency (rad/s) and sort
wavelength_m = wavelength_nm * 1e-9;
omega_rad_s = (2 * pi * c) ./ wavelength_m;
[omega, sortIdx] = sort(omega_rad_s);
I_omega = I_spectrum(sortIdx);

% --- APPLY TRANSFORM-LIMITED ASSUMPTION ---
% 1. Get Amplitude: Magnitude is the square root of the intensity.
Amplitude_omega = sqrt(I_omega);

% 2. Set Phase: For TL, the phase is zero.
Phase_omega = zeros(size(Amplitude_omega)); 

% 3. Create Complex Electric Field Spectrum E(omega)
E_omega_TL = Amplitude_omega .* exp(1i * Phase_omega); % Since exp(i*0)=1, this is just Amplitude_omega
N_original = length(E_omega_TL);

% --- 2. Zero Padding and Time Vector Setup ---

% Parameters for Time Domain Control (Define total size N)
T_total = 100*1e-15; % Total duration of the reconstructed pulse (e.g., 500 fs)
N_time = 2^14;     % Choose N to be a power of 2 for conceptual clarity (e.g., 4096)
dt = T_total / N_time;
t = (0:N_time-1)' * dt; % Time vector (N_time x 1)

% Perform Zero Padding
N_zeros = N_time - N_original;
if N_zeros < 0
    % Truncate and warn if N_time is too small
    E_padded = E_omega_TL(1:N_time); 
    N_eff = N_time;
else
    % Zero padding at the end
    E_padded = [E_omega_TL; zeros(N_zeros, 1)];
    N_eff = N_time;
end


% --- 3. Manual IDFT Calculation ---
E_t = zeros(N_eff, 1); % Initialize the resulting time-domain electric field

% Loop through each point in the time domain (n)
for n = 0 : N_eff - 1
    sum_val = 0;
    
    % Loop through each point in the frequency domain (k)
    for k = 0 : N_eff - 1
        % This is the Inverse Discrete Fourier Transform kernel
        exponent = 1i * (2*pi/N_eff) * k * n;
        sum_val = sum_val + E_padded(k + 1) * exp(exponent);
    end
    
    % Store the result for this time point
    % For a standard IDFT, the final scaling factor (1/N or similar) 
    % depends on the definition. We can apply it later for plotting.
    E_t(n + 1) = sum_val; 
end

% --- 4. Plotting and Analysis ---

% Shift the time axis for the pulse to be centered (optional, but typical)
E_t_shifted = fftshift(E_t);
t_shifted = t - T_total/2; 

I_t = abs(E_t_shifted).^2;

figure;
plot(t_shifted * 1e15, I_t / max(I_t)); % Plot time in femtoseconds (fs), normalized
xlabel('Time (fs)');
ylabel('Normalized Intensity (Arb. Units)');
title('Reconstructed Transform-Limited Pulse');