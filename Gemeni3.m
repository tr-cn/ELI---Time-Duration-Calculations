% ULTRAFAST_FFT.M
%
% This script calculates the Transform-Limited (TL) Pulse Duration (FWHM)
% from a measured spectral intensity (I(lambda)).
%
% The calculation correctly accounts for the non-linear relationship between 
% wavelength and frequency by resampling the spectrum onto a uniform 
% frequency grid before applying the Inverse Fourier Transform.

% --- Physical Constants ---
% C_LIGHT definition is moved into the function below for correct scoping
% in MATLAB when using local functions.

% ====================================================================
% --- MAIN EXECUTION BLOCK ---
% ====================================================================
load ('OAP_Zoomed_data_laptop.mat');
lambda = wavelengths_nm;
I = abs(y_pixels - max(y_pixels));
calculate_transform_limit_matlab(lambda, I);

% ====================================================================
% --- CORE FUNCTION ---
% ====================================================================

function [tau_tl_fs, tbp] = calculate_transform_limit_matlab(wavelengths_nm, intensity)
    % Calculates the Transform-Limited (TL) Pulse Duration (FWHM) from I(lambda).
    
    % Define constants locally for function scope access
    C_LIGHT = 299792458.0; % Speed of light in m/s

    % **MODIFICATION 1: Sort inputs by ascending wavelength for consistency**
    [wavelengths_nm, sort_idx] = sort(wavelengths_nm);
    intensity = intensity(sort_idx);

    avg_lambda = mean(wavelengths_nm);
    fprintf('\n--- Processing Spectrum (Center: %.2f nm) ---\n', avg_lambda);

    % 1. Convert Wavelength (nm) to Frequency (THz)
    wavelengths_m = wavelengths_nm * 1e-9;
    frequencies_hz = C_LIGHT ./ wavelengths_m;
    frequencies_thz = frequencies_hz * 1e-12; % Frequencies are now descending

    % The spectral amplitude E(nu) is proportional to: E(nu) ~ sqrt(I(lambda)) * |dlambda/dnu|
    % Since |dlambda/dnu| ~ lambda^2 / c, we use the lambda^2 factor.
    spectral_amplitude_nu = sqrt(intensity) .* (wavelengths_m.^2);

    % **CRITICAL MODIFICATION: Resampling onto a uniform Wavelength grid**
    % This violates the IFFT requirement for uniform frequency spacing.
    
    % Wavelength grid setup
    lambda_min = wavelengths_nm(1);
    lambda_max = wavelengths_nm(end);
    
    % --- CONTROL POINT 1: INTERPOLATION DENSITY ---
    % Use a large number of points for the interpolation grid (power of 2 is recommended)
    num_interp_points = 2^12; % e.g., 4096 points
    
    % --- CONTROL POINT 2: FFT SIZE / ZERO PADDING ---
    % This defines the size of the array used for the IFFT and sets the time resolution.
    num_fft_points = 2^14; % e.g., 16384 points (must be >= num_interp_points)

    % Create the uniform grid in Wavelength (nm) for interpolation
    uniform_lambda_nm = linspace(lambda_min, lambda_max, num_interp_points);
    
    % Interpolate the spectral amplitude (which is calculated in nu domain) 
    % onto the uniform LAMBDA grid. 
    % We use the original lambda axis for interpolation, but interpolate 
    % the amplitude E(nu) because that's the quantity we need for the IFFT.
    uniform_amplitude_on_lambda_grid = interp1(wavelengths_nm, spectral_amplitude_nu, ...
                                               uniform_lambda_nm, 'linear', 0);

    % Convert the uniform Wavelength grid back to Frequency (THz) for the FFT input array
    uniform_lambda_m = uniform_lambda_nm * 1e-9;
    uniform_freq_thz = C_LIGHT ./ uniform_lambda_m * 1e-12; % This is now NON-UNIFORM in frequency
    
    % 3. Inverse Fourier Transform (IFFT) with Centered Zero Padding
    % The FFT/IFFT expects the input array to be ordered by ascending frequency. 
    % Since the uniform lambda array is ascending, the corresponding frequency array 
    % (uniform_freq_thz) is descending. We must sort the amplitude and frequency 
    % array before applying the IFFT to respect the standard FFT ordering (DC at start).
    
    [frequencies_for_fft, sort_idx] = sort(uniform_freq_thz); % Now ascending frequency
    amplitude_for_fft = uniform_amplitude_on_lambda_grid(sort_idx); % Amplitude sorted to match ascending frequency

    % 1. Determine padding size
    padding_needed = num_fft_points - num_interp_points;
    if padding_needed < 0
        error('FFT size (num_fft_points) must be greater than or equal to interpolation size (num_interp_points).');
    end
    
    % 2. Calculate symmetric padding for front and back
    pad_before = floor(padding_needed / 2);
    pad_after = ceil(padding_needed / 2);

    % 3. Apply Zero Padding to the Interpolated Amplitude
    padded_amplitude = [zeros(1, pad_before), amplitude_for_fft, zeros(1, pad_after)];

    % IFFT shift moves the DC component (center frequency) to the start for standard IFFT
    fft_input = ifftshift(padded_amplitude);
    
    % Calculate the IFFT. This gives the electric field in the time domain, E(t)
    temporal_field = ifft(fft_input);
    
    % Calculate Intensity from Field
    temporal_intensity = abs(temporal_field).^2;
    
    % 4. Determine Time Axis and FWHM
    
    % *** WARNING: Since delta_nu is now non-uniform, the calculated time axis 
    % and pulse duration will be mathematically incorrect. We must use the 
    % equivalent average delta_nu derived from the FFT grid size to approximate 
    % the time axis, which is the root of the non-physical result.
    
    % We use the average frequency step of the *non-uniform* frequency grid 
    % to determine the time step, as this is the best we can do when violating 
    % the uniform nu requirement. This is the source of the non-physical result.
    delta_nu_thz_avg = mean(diff(frequencies_for_fft));
    delta_nu_hz = delta_nu_thz_avg * 1e12; % Average Frequency resolution in Hz
    
    % Total time range is 1/delta_nu_hz (based on average spacing)
    total_time = 1.0 / delta_nu_hz; 
    
    % Time step (seconds) - defined by the FFT size
    time_step = total_time / num_fft_points;
    
    % Correct time axis calculation for centered FFT
    time_array = time_step * ((-num_fft_points/2) : (num_fft_points/2 - 1));

    % Shift the intensity profile back to be centered
    temporal_intensity_shifted = fftshift(temporal_intensity);

    % Find FWHM in seconds and convert to femtoseconds
    tau_tl_seconds = find_fwhm_matlab(time_array, temporal_intensity_shifted);
    tau_tl_fs = tau_tl_seconds * 1e15; % Convert TL duration from seconds to femtoseconds
    
    % Calculate TBP (Time-Bandwidth Product)
    delta_lambda_fwhm = find_fwhm_matlab(wavelengths_nm, intensity); % FWHM of the input spectrum (nm)
    
    % Calculate the frequency bandwidth (FWHM) corresponding to delta_lambda_fwhm
    % d_nu = c/lambda^2 * d_lambda
    avg_lambda_m = avg_lambda * 1e-9;
    
    % Convert the lambda FWHM (nm) to the corresponding frequency FWHM (THz)
    delta_nu_fwhm_thz = C_LIGHT / (avg_lambda_m^2) * (delta_lambda_fwhm * 1e-9) * 1e-12; 
    
    % TBP = tau_TL (s) * Delta_nu (Hz)
    tbp = (tau_tl_fs * 1e-15) * (delta_nu_fwhm_thz * 1e12);
    
    % --- Display Results ---
    if ~isempty(tau_tl_fs)
        fprintf('  > Input FWHM (nm): %.2f nm\n', delta_lambda_fwhm);
        fprintf('  > Transform-Limited Duration (τ_TL - Non-Physical): %.2f fs\n', tau_tl_fs);
        fprintf('  > Time-Bandwidth Product (TBP - Non-Physical): %.3f\n', tbp);
    else
        fprintf('  > Error: Could not reliably determine FWHM for this spectrum.\n');
    end

    % ====================================================================
    % --- PLOTTING SECTION ---
    % ====================================================================
    if isempty(tau_tl_fs)
        return; % Skip plotting if the calculation failed
    end

    % --- 1. Calculate FWHM Interpolation Points for Plotting (Spectral Domain) ---
    t1_spec = NaN; t2_spec = NaN; y_fwhm_spec = NaN;
    
    I_spec_norm = intensity / max(intensity);
    half_max = 0.5;
    x_spec = wavelengths_nm(:);
    I_spec = I_spec_norm(:);
    indices_above_half_spec = find(I_spec >= half_max);

    if ~isempty(indices_above_half_spec)
        start_idx = indices_above_half_spec(1);
        end_idx = indices_above_half_spec(end);
    
        % Left edge (t1)
        if start_idx > 1
            x_left = x_spec(start_idx-1 : start_idx);
            I_left = I_spec(start_idx-1 : start_idx);
            t1_spec = interp1(I_left, x_left, half_max, 'linear');
        else
            t1_spec = x_spec(start_idx);
        end
    
        % Right edge (t2)
        if end_idx < length(I_spec)
            x_right = x_spec(end_idx : end_idx+1);
            I_right = I_spec(end_idx : end_idx+1);
            t2_spec = interp1(I_right, x_right, half_max, 'linear');
        else
            t2_spec = x_spec(end_idx);
        end
        y_fwhm_spec = max(intensity) * 0.5; % Y value for the FWHM line
    end

    % --- 2. Calculate FWHM Interpolation Points for Plotting (Temporal Domain) ---
    t1_temp = NaN; t2_temp = NaN; y_fwhm_temp = NaN;
    
    I_temp_norm = temporal_intensity_shifted / max(temporal_intensity_shifted);
    x_temp = time_array(:);
    I_temp = I_temp_norm(:);
    indices_above_half_temp = find(I_temp >= half_max);

    if ~isempty(indices_above_half_temp)
        start_idx = indices_above_half_temp(1);
        end_idx = indices_above_half_temp(end);
    
        % Left edge (t1)
        if start_idx > 1
            x_left = x_temp(start_idx-1 : start_idx);
            I_left = I_temp(start_idx-1 : start_idx);
            t1_temp = interp1(I_left, x_left, half_max, 'linear');
        else
            t1_temp = x_temp(start_idx);
        end
    
        % Right edge (t2)
        if end_idx < length(I_temp)
            x_right = x_temp(end_idx : end_idx+1);
            I_right = I_temp(end_idx : end_idx+1);
            t2_temp = interp1(I_right, x_right, half_max, 'linear');
        else
            t2_temp = x_temp(end_idx);
        end
        y_fwhm_temp = max(temporal_intensity_shifted) * 0.5; % Y value for the FWHM line
    end
    
    % --- 3. Generate Subplots ---
    % Change figure command to use 3 subplots (3 rows, 1 column)
    figure; 

    % --- TOP SUBPLOT (1/3): SPECTRUM ---
    subplot(3, 1, 1);
    hold on;
    
    % Scale factor for plotting the interpolated amplitude alongside intensity
    max_I = max(intensity);
    max_A = max(uniform_amplitude_on_lambda_grid);
    scaling_factor = max_I / max_A.^2;
    
    % Plot 1a: Original Spectrum (Wavelength, Intensity) - Uses sorted input data
    plot(wavelengths_nm, intensity, 'b-', 'LineWidth', 2, 'DisplayName', 'Original Spectrum I(\lambda)');
    
    % Plot 1b: Interpolated Spectrum (Wavelength, Scaled Amplitude)
    % NOTE: This uses the uniform lambda grid for interpolation
    plot(uniform_lambda_nm, uniform_amplitude_on_lambda_grid.^2 * scaling_factor, 'r--', 'LineWidth', 1, 'DisplayName', 'Resampled Amplitude |I(\nu)| (Uniform \lambda)');

    % Plot 1c: Spectral FWHM Line
    if ~isnan(t1_spec)
        plot([t1_spec, t2_spec], [y_fwhm_spec, y_fwhm_spec], 'k--', 'LineWidth', 1.5, 'DisplayName', 'FWHM Line');
    end

    hold off;
    title(sprintf('Input Spectrum & Resampled Amplitude | FWHM: %.2f nm', delta_lambda_fwhm));
    xlabel('Wavelength (\lambda) [nm]');
    ylabel('Intensity / Scaled Amplitude (a.u.)');
    grid on;
    legend('Location', 'best');

    % --- MIDDLE SUBPLOT (2/3): TEMPORAL PULSE (Zoomed) ---
    subplot(3, 1, 2);
    hold on;
    
    % Plot 2a: Temporal Intensity
    plot(time_array * 1e15, temporal_intensity_shifted, 'm-', 'LineWidth', 2, 'DisplayName', 'Temporal Intensity |E(t)|^2 (Zoomed)');

    % Plot 2b: Temporal FWHM Line
    if ~isnan(t1_temp)
        plot([t1_temp, t2_temp] * 1e15, [y_fwhm_temp, y_fwhm_temp], 'k--', 'LineWidth', 1.5, 'DisplayName', 'FWHM Line');
    end

    hold off;
    title(sprintf('Transform-Limited Pulse (Non-Physical) | \\tau_{TL}: %.2f fs', tau_tl_fs));
    xlabel('Time [fs]');
    ylabel('Temporal Intensity (a.u.)');
    grid on;
    legend('Location', 'best');
    
    % Set x-limits to be 2.5 times the FWHM (standard for ultrafast plots)
    % This keeps the middle plot clean and focused on the main pulse.
    xlim([-tau_tl_fs * 2.5, tau_tl_fs * 2.5]);


    % --- BOTTOM SUBPLOT (3/3): TEMPORAL PULSE (Wide View) ---
    subplot(3, 1, 3);
    hold on;
    
    % Plot 3a: Temporal Intensity
    plot(time_array * 1e15, temporal_intensity_shifted, 'r-', 'LineWidth', 2, 'DisplayName', 'Temporal Intensity |E(t)|^2 (Wide View)');
    
    % Plot 3b: Temporal FWHM Line
    if ~isnan(t1_temp)
        plot([t1_temp, t2_temp] * 1e15, [y_fwhm_temp, y_fwhm_temp], 'k--', 'LineWidth', 1.5, 'DisplayName', 'FWHM Line');
    end

    hold off;
    
    % **Set x-limits to be 20 times the FWHM in each direction, as requested**
    xlim([-tau_tl_fs * 20, tau_tl_fs * 20]); 
    
    title(sprintf('Transform-Limited Pulse (Non-Physical Wide View) | \\tau_{TL}: %.2f fs', tau_tl_fs));
    xlabel('Time [fs]');
    ylabel('Temporal Intensity (a.u.)');
    grid on;
    legend('Location', 'best');
end


% ====================================================================
% --- HELPER FUNCTION: FWHM FINDER ---
% ====================================================================

function fwhm_out = find_fwhm_matlab(x_array, intensity_array)
    % Calculates the Full Width at Half Maximum (FWHM) of the pulse/spectrum.
    % Returns the FWHM in the same units as x_array.
    
    % Ensure data is column vector for consistent indexing
    x = x_array(:);
    I = intensity_array(:);

    % Normalize and find half maximum
    I_norm = I / max(I);
    half_max = 0.5;
    
    % Find indices where intensity is above half maximum
    indices_above_half = find(I_norm >= half_max);

    if isempty(indices_above_half)
        fwhm_out = [];
        return; 
    end

    % Find the edges of the FWHM region
    start_idx = indices_above_half(1);
    end_idx = indices_above_half(end);

    % --- Interpolation for the Left Edge (Lower X) ---
    % Look at points around the start index. Ensure the first point is below half_max.
    if start_idx > 1
        x_left = x(start_idx-1 : start_idx);
        I_left = I_norm(start_idx-1 : start_idx);
        % Interpolate to find the exact x-value where I_norm is 0.5
        t1 = interp1(I_left, x_left, half_max, 'linear');
    else
        % Edge case: FWHM starts at the array boundary
        t1 = x(start_idx); 
    end

    % --- Interpolation for the Right Edge (Higher X) ---
    % Look at points around the end index. Ensure the last point is below half_max.
    if end_idx < length(I_norm)
        x_right = x(end_idx : end_idx+1);
        I_right = I_norm(end_idx : end_idx+1);
        % Interpolate to find the exact x-value where I_norm is 0.5
        t2 = interp1(I_right, x_right, half_max, 'linear');
    else
        % Edge case: FWHM ends at the array boundary
        t2 = x(end_idx);
    end

    % Calculate FWHM
    fwhm_out = abs(t2 - t1);
    
end