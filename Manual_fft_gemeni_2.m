clear; close all; clc;
%% Inverse Fourier Transform for Pulse Reconstruction (Corrected for 2*pi scaling)
% --- 1. Load and Convert Data ---
% load("OAP_Zoomed_data_laptop.mat")
% load("data.mat")
load("Real_data.mat")
c = 299792458; % Speed of light (m/s)

% Your Input Data (Wavelength in nm, Spectrum is REAL Intensity)
wavelength_nm = wavelengths_nm; % Example: Your input vector in nm
I_spectrum_meas = abs(y_pixels - max(y_pixels));
[wavelength_nm,indx]  = sort(wavelength_nm);
I_spectrum_meas = I_spectrum_meas(indx);

% --- 2. USER INPUT: DESIRED TIME DOMAIN CONTROL ---
N_eff = 2^20;        % The desired number of points for the IFFT (controls T_total length)

fprintf('\n--- Input Parameters ---\n');
fprintf('FFT Points (N): %d\n', N_eff);

% --- 3. CONVERSION AND RESAMPLING TO UNIFORM FREQUENCY GRID ---

% 3.1. Convert Wavelength to non-uniform Frequency (Hz)
wavelength_m = wavelength_nm * 1e-9;
nu_Hz = c ./ wavelength_m; 
E_nu_TL_non_uniform = sqrt(I_spectrum_meas); % Transform-Limited E-field (Amplitude)

% 3.2. Define the new, UNIFORM Frequency Grid
% Sort in ascending frequency order
[nu_Hz, sortIdx] = sort(nu_Hz);
E_nu_TL_non_uniform = E_nu_TL_non_uniform(sortIdx);

% Determine the step size for the uniform grid (use the original number of points over the span)
nu_span = nu_Hz(end) - nu_Hz(1);
dnu_uniform = nu_span / (length(nu_Hz) - 1); 

% Create the new uniform grid spanning the original data range
nu_uniform_grid = linspace(nu_Hz(1),nu_Hz(end),2^10);%nu_Hz(1) : dnu_uniform : nu_Hz(end);

% 3.3. Resample the spectrum onto the uniform grid
% Use 'linear' interpolation, setting values outside the range to zero (0).
% Transpose to column vector (required for subsequent operations)
E_nu_resampled = interp1(nu_Hz, E_nu_TL_non_uniform, nu_uniform_grid, 'linear', 0).'; 

% --- 4. PARAMETER CALCULATION (Driven by Uniform Grid) ---

% 4.1. Calculate the IFFT Total Time Range (T_total)
% T_total is dictated by the uniform frequency step: T_total = 1 / Delta_nu_eff
Delta_nu_eff = nu_uniform_grid(2) - nu_uniform_grid(1);
T_total = 1 / Delta_nu_eff; 

% 4.2. Calculate the resulting Time Step (dt)
dt_actual = T_total / N_eff; 

fprintf('Resampled Frequency Step (Delta nu): %.2f GHz\n', Delta_nu_eff / 1e9);
fprintf('Total Time Range (T_total): %.2f fs\n', T_total * 1e15);
fprintf('Calculated Time Step (dt): %.2e fs\n', dt_actual * 1e15);

% --- 5. CENTERING AND PADDING ---

N_resampled = length(E_nu_resampled);

% 5.1. Center the spectrum (moves carrier to the center of the vector)
% E_shifted = ifftshift(E_nu_resampled);

% 5.2. Apply Padding and Un-shift
N_shift = floor((N_eff - N_resampled) / 2); 
E_padded_centered = [zeros(N_shift, 1); E_nu_resampled; zeros(N_eff - N_resampled - N_shift, 1)];
E_IFFT_input = ifftshift(E_padded_centered); % Un-shift for IFFT (low frequency at index 1)

% --- 6. EXECUTE INVERSE FAST FOURIER TRANSFORM (IFFT) ---
E_t = fft(E_IFFT_input); 
E_t = fftshift(E_t);
% --- 7. TIME VECTOR CONSTRUCTION AND ANALYSIS ---

t_raw = (0:N_eff-1)' * dt_actual;
t_shifted = t_raw - T_total/2; % Center the time axis around zero

I_t = abs(E_t).^2;
I_t_norm = I_t / max(I_t);

% Calculate FWHM
[I_peak, idx_peak] = max(I_t_norm);
half_max = 0.5 * I_peak; % Use 0.5 * actual peak value

% --- Robust FWHM Calculation using Two-Point Linear Interpolation ---
if I_peak < half_max
    warning('Pulse peak is below half-maximum. FWHM calculation aborted.');
    FWHM = NaN; t1 = NaN; t2 = NaN;
else
    % --- Find t1 (Leading Edge) ---
    I_before = I_t_norm(1:idx_peak);
    t_before = t_shifted(1:idx_peak);
    
    idx_start_t1 = find(I_before <= half_max, 1, 'last'); % Point <= 0.5
    idx_end_t1 = find(I_before >= half_max, 1, 'first');  % Point >= 0.5
    
    if isempty(idx_start_t1) || isempty(idx_end_t1)
        warning('Cannot reliably find leading edge (t1) points.');
        t1 = NaN;
    else
        % Use the two surrounding points for linear interpolation
        t1_indices = [idx_start_t1, idx_end_t1];
        I_segment_t1 = I_before(t1_indices);
        t_segment_t1 = t_before(t1_indices);
        
        [I_unique_t1, unique_idx_t1] = unique(I_segment_t1);
        t1 = interp1(I_unique_t1, t_segment_t1(unique_idx_t1), half_max, 'linear');
    end

    % --- Find t2 (Trailing Edge) ---
    I_after = I_t_norm(idx_peak:end);
    t_after = t_shifted(idx_peak:end);

    idx_start_t2 = find(I_after >= half_max, 1, 'last');  % Point >= 0.5
    idx_end_t2 = find(I_after <= half_max, 1, 'first');   % Point <= 0.5

    if isempty(idx_start_t2) || isempty(idx_end_t2)
        warning('Cannot reliably find trailing edge (t2) points.');
        t2 = NaN;
    else
        % Use the two surrounding points for linear interpolation
        t2_indices = [idx_start_t2, idx_end_t2];
        I_segment_t2 = I_after(t2_indices);
        t_segment_t2 = t_after(t2_indices);

        [I_unique_t2, unique_idx_t2] = unique(I_segment_t2);
        t2 = interp1(I_unique_t2, t_segment_t2(unique_idx_t2), half_max, 'linear');
    end
    
    FWHM = t2 - t1;
end 

fprintf('Calculated FWHM Duration: %.2f fs\n', FWHM * 1e15);

% --- 8. PLOTTING ---
figure;
subplot(2,1,1)
plot(wavelength_nm,I_spectrum_meas, "LineWidth",2,"DisplayName","Measured Data"); hold on
plot(c./nu_uniform_grid*1e9,abs(E_nu_resampled).^2,"LineWidth",2,"LineStyle","--","Color", "red","DisplayName","Intepulated Data")
legend


subplot(2,1,2)
plot(t_shifted * 1e15, I_t_norm, 'LineWidth', 2);
hold on;
% Plot the FWHM line
if ~isnan(FWHM)
    plot([t1 t2] * 1e15, [half_max half_max], 'r--', 'LineWidth', 1.5);
    plot([t1 t1] * 1e15, [0 half_max], 'r:');
    plot([t2 t2] * 1e15, [0 half_max], 'r:');
end
hold off;

xlabel('Time (fs)', 'FontSize', 14);
ylabel('Normalized Intensity', 'FontSize', 14);
title(['Reconstructed TL Pulse (FWHM: ' num2str(FWHM*1e15, '%.2f') ' fs)'], 'FontSize', 16);
grid on;
axis([min(t_shifted)*1e15 max(t_shifted)*1e15 0 1.05]); % Set y-axis limit
xlim([-150,150])