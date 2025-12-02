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
    points_x = [];
    points_y = [];
    setappdata(gcf, 'points_x', points_x);
    setappdata(gcf, 'points_y', points_y);
    setappdata(gcf, 'digitizing', true);
    
    % Step 3: Set up callbacks
    set(gcf, 'WindowButtonDownFcn', @onClick);
    set(gcf, 'KeyPressFcn', @onKeyPress);
    
    % Instructions
    text(0.02, 0.98, {'Left-click to digitize points', ...
                      'Press ESC to finish & analyze'}, 'Units', 'normalized', ...
                      'VerticalAlignment', 'top', 'FontSize', 12);
    
    % Step 4: Wait for user to finish
    uiwait(gcf);
    
    % Step 5: Get digitized data
    digitized_x = getappdata(gcf, 'points_x');
    digitized_y = getappdata(gcf, 'points_y');
    
    % Step 6: Save to workspace
    assignin('base', 'digitized_x', digitized_x);
    assignin('base', 'digitized_y', digitized_y);
    
    % Step 7: Analyze with Fourier transform (using our custom functions)
    fprintf('\n=== Digitized %d points ===\n', length(digitized_x));
    
    % Forward DFT: time -> frequency
    Ew = myDFT(digitized_y);
    N = length(digitized_y);
    dt = digitized_x(2) - digitized_x(1);  % time step
    domega = 2*pi/(N*dt);                  % frequency spacing
    
    % Frequency axis
    omega = (0:N-1)' * domega;
    
    % Plot frequency domain
    figure('Position', [200, 200, 1000, 400]);
    subplot(1,2,1);
    plot(omega/(2*pi), abs(Ew)/N);  % magnitude spectrum
    title('Frequency Spectrum |E(ω)|'); xlabel('Frequency (Hz)'); grid on;
    
    % Inverse DFT: frequency -> time (verify it works)
    y_recovered = myIDFT(Ew);
    
    subplot(1,2,2);
    plot(digitized_x, digitized_y, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    hold on;
    plot(digitized_x, real(y_recovered), 'b-', 'LineWidth', 2);
    legend('Digitized', 'Recovered (IDFT)', 'Location', 'best');
    title('Original vs Recovered'); xlabel('Time'); ylabel('Amplitude'); grid on;
    
    % Error stats
    error = max(abs(digitized_y - real(y_recovered)));
    fprintf('Reconstruction error: %.2e (should be ~machine precision)\n', error);
    
    fprintf('Data saved to workspace: digitized_x, digitized_y\n');
end

% === OUR CUSTOM DFT FUNCTIONS ===
function X = myDFT(x)
    x = x(:); N = length(x);
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
