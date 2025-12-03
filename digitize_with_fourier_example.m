function digitize_with_fourier_example
    % Step 1: Create or load a figure to digitize (example noisy sine wave)
    % figure('Position', [100, 100, 800, 600]);
    % t_original = 0:0.1:10;
    % y_original = sin(2*pi*t_original) + 0.1*randn(size(t_original));  % noisy
    % plot(t_original, y_original, 'b-', 'LineWidth', 2); hold on; grid on;
    % title('Digitize this curve (ESC to stop)');
    
    % figure('Position', [100, 100, 800, 600]);
    % img = imread('Real spectrum.jpeg');
    % imshow(img)



    % Step 2: Initialize digitizer
    % points_x = [];
    % points_y = [];
    % setappdata(gcf, 'points_x', points_x);
    % setappdata(gcf, 'points_y', points_y);
    % setappdata(gcf, 'digitizing', true);
    % 
    % % Step 3: Set up callbacks
    % set(gcf, 'WindowButtonDownFcn', @onClick);
    % set(gcf, 'KeyPressFcn', @onKeyPress);
    % 
    % % Instructions
    % text(0.02, 0.98, {'Left-click to digitize points', ...
    %                   'Press ESC to finish & analyze'}, 'Units', 'normalized', ...
    %                   'VerticalAlignment', 'top', 'FontSize', 12);
    % 
    % % Step 4: Wait for user to finish
    % uiwait(gcf);
    % 
    % % Step 5: Get digitized data
    % digitized_x = getappdata(gcf, 'points_x');
    % digitized_y = getappdata(gcf, 'points_y');
    % 
    % % Step 6: Save to workspace
    % assignin('base', 'digitized_x', digitized_x);
    % assignin('base', 'digitized_y', digitized_y);
    % 
    % % Step 7: Analyze with Fourier transform (using our custom functions)
    % fprintf('\n=== Digitized %d points ===\n', length(digitized_x));
    % 
    % % Forward DFT: time -> frequency
    load ("Real_data.mat")
    N_interp = 2^10;
    padding = 2;
    wavelengths_lin_nm = linspace(wavelengths_nm(1),wavelengths_nm(end),N_interp);
    dlabmda_nm = wavelengths_lin_nm(2) - wavelengths_lin_nm(1);
    
    y_pixels = abs(y_pixels-max(y_pixels));
    
    y_pixels_interp = interp1(wavelengths_nm,y_pixels,wavelengths_lin_nm);
    
    Et = myDFT(y_pixels_interp,padding);
    Et_centerd = fftshift(Et);
    N = length(y_pixels_interp);
    % dlabmda_nm = (wavelengths_nm(2) - wavelengths_nm(1))/N;  % time step
    
    
    dF_hz = 3e8/(dlabmda_nm*1e-9);
    dt_s = 1/(N*dF_hz);                  % time spacing
    
    % Frequency axis
    t_s = (0:(N+N*padding)-1)' * dt_s;
    % wavelengths_lin_padded_nm = linspace(wavelengths_lin_nm(1),wavelengths_lin_nm(end),length(t_s));
    wavelengths_lin_padded_nm = [zeros(1,N*padding/2), wavelengths_lin_nm,zeros(1,N*padding/2)];
    % Plot frequency domain
    figure('Position', [200, 200, 1000, 400]);
    subplot(1,2,1);
    plot(t_s/1e-15, abs(Et_centerd)/N);  % magnitude spectrum
    title('Frequency Spectrum |E(t)|'); xlabel('time (fs)'); grid on;
    
    % Inverse DFT: frequency -> time (verify it works)
    y_recovered = myIDFT(Et);
    
    subplot(1,2,2);
    plot(wavelengths_lin_nm, y_pixels_interp, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    hold on;
    plot(wavelengths_lin_padded_nm, real(y_recovered), 'b-', 'LineWidth', 2);
    legend('Digitized', 'Recovered (IDFT)', 'Location', 'best');
    title('Original vs Recovered'); xlabel('Wavelength (nm)'); ylabel('Amplitude'); grid on;
    
    % Error stats
    error = max(abs(y_pixels_interp - real(y_recovered)));
    % fprintf('Reconstruction error: %.2e (should be ~machine precision)\n', error);
    
    % fprintf('Data saved to workspace: wavelengths_nm, y_pixels\n');
end

% === OUR CUSTOM DFT FUNCTIONS ===
function X = myDFT(x,padding)
    N = length(x) + length(x)*padding ; x = [zeros(length(x)*padding/2,1); x(:); zeros(length(x)*padding/2,1)];
    X = zeros(N,1);
    for k = 0:N-1
        s = 0;
        for n = 0:N-1
            s = s + x(n+1)*exp(-1j*2*pi*k*n/N);
        end
        X(k+1) = s;
    end
end

function e = myIDFT(Ew)
    Ew = Ew(:); N = length(Ew);
    e = zeros(N,1);
    for n = 0:N-1
        s = 0;
        for k = 0:N-1
            s = s + Ew(k+1)*exp(1j*2*pi*k*n/N);
        end
        e(n+1) = s / N;
    end
end

% === DIGITIZER CALLBACKS ===
function onClick(~, ~)
    if ~getappdata(gcf, 'digitizing'), return; end
    
    cp = get(gca, 'CurrentPoint');
    x = cp(1,1); y = cp(1,2);
    
    points_x = getappdata(gcf, 'points_x');
    points_y = getappdata(gcf, 'points_y');
    points_x(end+1) = x;
    points_y(end+1) = y;
    setappdata(gcf, 'points_x', points_x);
    setappdata(gcf, 'points_y', points_y);
    
    % Show dot immediately
    plot(x, y, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'LineWidth', 2);
    title(sprintf('Digitized (%d pts), ESC to analyze with Fourier', length(points_x)));
    drawnow;
end

function onKeyPress(~, event)
    if strcmp(event.Key, 'escape')
        setappdata(gcf, 'digitizing', false);
        uiresume(gcf);
    end
end
