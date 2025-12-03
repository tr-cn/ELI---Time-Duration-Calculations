% clear; close all; clc;
if false
    % clear; close all; clc;
    % img = imread('OAP spectrum Zoomed.jpeg');
    % img = imread('Real spectrum.jpeg');
    
    figure;
    imshow(img);
    imshow(imrotate(img,0.3));% OAP +0.3 % Real -1
    point = 40;
    title(['Click ',num2str(point),' points'])
    [x_pixels, y_pixels] = ginput(point);
    
    
    hold on
    scatter(x_pixels,y_pixels,'filled','b')


    x_ref1_pixel = 83;%OAP Zoomed % OAP 278;% Real 220;%
    x_ref1_value = 750;% 650;%

    x_ref2_pixel = 1230;%OAP_Real%OAP855;% 1145;%
    x_ref2_value = 900;% 1000;%

    m_wavelength = (x_ref2_value - x_ref1_value) / (x_ref2_pixel - x_ref1_pixel);
    b_wavelength = x_ref1_value - m_wavelength * x_ref1_pixel;

    wavelengths_nm = m_wavelength * x_pixels + b_wavelength;
    disp('--- Digitized Data ---');
    disp(['Wavelength (nm): ', num2str(wavelengths_nm')]);
    disp(['Y-Values (pixels): ', num2str(y_pixels')]); % ‘‚Ë€Ÿ› ‘·’‰ŸŸ› Ÿ‘Ÿ’ ‹–◊Ë €Ÿ’‹ Y

    save("data.mat")
end

if true
% clear; close all; clc
load("Real_data.mat")
% load("Real_data.mat")
c = 299792458; % m
interp_steps = 2^12;
FFT_steps = 2^20;

[wavelengths_nm, indx] = sort(wavelengths_nm);
wavelengths_m = wavelengths_nm*1e-9;
y_pixels = y_pixels(indx);

I_m = abs(y_pixels - max(y_pixels)); % turning the - sign to pulse sign
E_m = sqrt(I_m);


freq_Hz = c./wavelengths_m;
freq_PHz = freq_Hz / 1e15;
omega_PHz = 2 * pi * freq_PHz;

omega_PHz_lin = linspace(omega_PHz(1),omega_PHz(end),interp_steps);
E_PHz    = interp1(omega_PHz, E_m, omega_PHz_lin, 'pchip', 0);
% E_PHz = E_PHz .* exp(pi*1*(randn(size(E_PHz))-0.0));

dt = 1/abs(freq_PHz(end) - freq_PHz(1)) / (FFT_steps/interp_steps); % fs

time_fs  = [-FFT_steps/2 : (FFT_steps/2-1) ] * dt;

E_fs_uncentered = ifft(E_PHz,FFT_steps);
figure; plot(real(E_fs_uncentered));
E_fs = ifftshift(E_fs_uncentered);
plot(time_fs,real(E_fs).^2);
I_fs = abs(E_fs).^2;
I_fs = I_fs/max(I_fs);
I_half_fs = 0.5;
indx_half = find(I_fs>=0.5);
fwhm_left_indx = indx_half(1);
fwhm_right_indx = indx_half(end);
fwhm_fs = abs(time_fs(fwhm_right_indx)-time_fs(fwhm_left_indx));

% disp(fwhm_fs)
disp('===================================================');
fprintf('Calculated Time-Limited Pulse Duration (FWHM): %.2f fs\n', fwhm_fs);
disp('===================================================');

figure;

subplot(2,1,1);
plot(wavelengths_nm, I_m, 'k', 'LineWidth', 2); hold on
% plot()
title('Extracted Spectral Intensity vs. Wavelength');
xlabel('Wavelength (nm)');
ylabel('Intensity (a.u.)');
grid on;

subplot(2,1,2);
plot(time_fs, I_fs, 'r', 'LineWidth', 2);
hold on;
plot([fwhm_left_indx fwhm_right_indx], [I_half_fs I_half_fs], 'b--', 'LineWidth', 1.5); % Á’ FWHM
title('Time-Limited Pulse Intensity Profile');
xlabel('Time (fs)');
ylabel('Normalized Intensity');
legend('Pulse Intensity', ['FWHM = ' num2str(fwhm_fs, '%.2f') ' fs'], 'Location', 'best');
xlim([-5*fwhm_fs, 5*fwhm_fs]); % ‘“—‹Í ‘ÊŸË ‹ÿ’’◊ Ë‹’’‡ÿŸ
grid on;




end



