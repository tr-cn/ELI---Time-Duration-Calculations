if false
    % щмб 1: итйръ дъоерд
    % дзмйфе аъ 'your_image.png' бщн дчебх щмлн
    img = imread('OAP spectrum.jpeg');
    % img = imread('Real spectrum.jpeg');
    
    figure;
    imshow(img);
    imshow(imrotate(img,+0.3));
    title('Click 25 points')

    % щмб 2: бзйшд йгрйъ щм рчегеъ одвшу баоцтеъ дтлбш
    % ginput оафщш бзйшъ рчегеъ вшфйъ
    [x_pixels, y_pixels] = ginput(30);
    hold on
    scatter(x_pixels,y_pixels,'filled','b')

    % x_pixels е-y_pixels олймйн лтъ аъ чеаешгйриеъ дфйчсмйн щм дрчегеъ щбзшъ.

    % щмб 3: лйем - дошд офйчсмйн мтшлйн фйжййн (аешк вм брреоишйн)
    % мцешк гевод, ррйз щдцйш даефчй деа аешк вм:

    % 1. двгшъ рчегеъ ййзес йгрйеъ (оъек цйш д-X щм двшу бъоерд)
    % йщ моцеа щъй рчегеъ йгетеъ тм дцйш даефчй, мгевод:
    x_ref1_pixel = 272; % дфйчсм щм ъзймъ дцйш даефчй
    x_ref1_value = 750; % дтшк дфйжй брчегд же (мгевод: 300 рреоиш)

    x_ref2_pixel = 856; % дфйчсм щм сеу дцйш даефчй
    x_ref2_value = 900; % дтшк дфйжй брчегд же (мгевод: 800 рреоиш)

    % 2. зйщеб фшоишй дишрсфешоцйд дмйрйашйъ:
    % ррйз чщш мйрйашй: Wavelength = m * X_pixel + b
    m_wavelength = (x_ref2_value - x_ref1_value) / (x_ref2_pixel - x_ref1_pixel);
    b_wavelength = x_ref1_value - m_wavelength * x_ref1_pixel;

    % 3. дошъ чеаешгйриеъ д-X щм двшу мтшлй аешк вм:
    wavelengths_nm = m_wavelength * x_pixels + b_wavelength;

    % (йщ мбцт ъдмйк геод тбеш цйш д-Y мдошъ тшлй дферчцйд)
    % ... (Y-axis calibration code here) ...

    % щмб 4: дцвъ дъецад
    disp('--- Digitized Data ---');
    disp(['Wavelength (nm): ', num2str(wavelengths_nm')]);
    disp(['Y-Values (pixels): ', num2str(y_pixels')]); % дтшлйн дсефййн йдйе мазш лйем Y

    % рйъп мщоеш аъ дръерйн мчебх:
    % dlmwrite('digitized_data.txt', [wavelengths_nm, y_values_calibrated]);
    save("data.mat")
end

if true
clear; close all; clc
load("OAP_data.mat")
% =======================================================
% --- щмб 1: двгшъ чбетйн еитйръ ръерйн (йщ мтглп!) ---
% =======================================================

% чбет одйшеъ даеш
% щйоещ бйзйгеъ резеъ маефийчд щм фемсйн чцшйн: рреоиш/фоиещрййд (nm/fs)
c = 299.792458; % c H 300 nm/fs 

% дзму аъ дгевоаеъ дмме бръерйн даойъййн щмк одгйвйицйд:
% (ега щаешлй двм осегшйн бсгш темд)
% -----------------------------------------------------------------
% 1. аешлй вм щргвое (одцйш даефчй, б-nm)
wavelengths_nm =wavelengths_nm; % <-- щрд аъ жд!

% [wavelengths_nm, indx] = sort(wavelengths_nm);
% 2. тецод сфчишмйъ щргвод (одцйш дарлй)
y_pixels =  abs(y_pixels - max(y_pixels));
Y_spectrum = y_pixels; % <-- щрд аъ жд!
% -----------------------------------------------------------------
% фшоишй гвйод м-FFT:
N_points = 2^17; % осфш рчегеъ FFT бсйсй (зжчд щм 2)
P_factor = 100;    % вешн Zero Padding (овгйм шжемецййъ жоп)
N_points_Padded = N_points * P_factor; % осфш рчегеъ FFT сефй

% 1. зйщеб дъгйшеъ джеейъйъ (omega) бйзйгеъ 1/fs
omega_orig = 2 * pi * c ./ wavelengths_nm; 

% 2. зйщеб ощштъ дщгд (|E(omega)|)
E_omega_magnitude_orig = sqrt(Y_spectrum);

% 3. ойеп: FFT гешщ сгш темд щм цйш дгвйод, млп ооййрйн аъ omega
[omega_orig_sorted, sort_idx] = sort(omega_orig);
E_omega_magnitude_orig_sorted = E_omega_magnitude_orig(sort_idx);




% 1. чбйтъ шщъ д-omega дзгщд (мйрйашйъ)
omega_min = min(omega_orig_sorted);
omega_max = max(omega_orig_sorted);
omega_lin = linspace(omega_min, omega_max, N_points);

% 2. айришфемцйд щм ощштъ дщгд мшщъ дмйрйашйъ
E_omega_TL_Interpolated = interp1(omega_orig_sorted, E_omega_magnitude_orig_sorted, omega_lin, 'pchip', 0);
% TL - Time-Limited (орйзйн фажд афс)

% 3. ййщен Zero Padding
% йецшйн ечиеш вгем йеъш (N_points_Padded) еоочойн аъ дръерйн доаеришфмйн бзмч дщоамй
E_omega_Padded = zeros(1, N_points_Padded);
E_omega_Padded(1:N_points) = E_omega_TL_Interpolated;




% 1. чбйтъ шщъ д-omega дзгщд (мйрйашйъ)
omega_min = min(omega_orig_sorted);
omega_max = max(omega_orig_sorted);
omega_lin = linspace(omega_min, omega_max, N_points);

% 2. айришфемцйд щм ощштъ дщгд мшщъ дмйрйашйъ
E_omega_TL_Interpolated = interp1(omega_orig_sorted, E_omega_magnitude_orig_sorted, omega_lin, 'pchip', 0);
% TL - Time-Limited (орйзйн фажд афс)

% йцйшъ змеп Hann баешк N_points
hann_window = hann(N_points)'; 

% длфмъ дсфчишен доаеришфм бзмеп
E_omega_windowed = E_omega_TL_Interpolated .* hann_window;

% 3. ййщен Zero Padding (ощъощйн бръерйн доезмрйн)
E_omega_Padded = zeros(1, N_points_Padded);
E_omega_Padded(1:N_points) = E_omega_windowed;

% 3. ййщен Zero Padding
% йецшйн ечиеш вгем йеъш (N_points_Padded) еоочойн аъ дръерйн доаеришфмйн бзмч дщоамй
E_omega_Padded = zeros(1, N_points_Padded);
E_omega_Padded(1:N_points) = E_omega_TL_Interpolated;



% 1. зйщеб цтг джоп (dt) догейч
% цтг джоп рчбт т"й иеез дъгшйн длемм евегм д-FFT (N_points_Padded).
Omega_Span = omega_max - omega_min; % [1/fs]
dt = (2 * pi / Omega_Span) / N_points_Padded * N_points; % [fs]

% цйш жоп ооешлж сбйб 0, бйзйгеъ фоиещрйеъ (fs)
     = (-N_points_Padded/2 : N_points_Padded/2 - 1) * dt;

% 2. бйцет IFFT
E_time_uncentered = ifft(E_omega_Padded);

E_time = fftshift(E_time_uncentered);

% 3. зйщеб тецоъ дфемс бжоп ешоем
I_time = abs(E_time).^2; 
I_time_norm = I_time / max(I_time);

% 4. оцйаъ ощк дфемс (FWHM)
I_half = 0.5;
% half_max_indices = find(I_time_norm >= I_half);

[~,indx]= max(I_time_norm);
half_max_indice_left = interp1(I_time_norm(1:indx),1:indx,0.5);
half_max_indice_right = interp1(I_time_norm(indx:end),indx:length(I_time_norm),0.5);

if isempty(half_max_indice_left || half_max_indice_right)
    tau_TL = NaN;
else
    % айргчсйн мдъзмд емсеу д-FWHM
    idx_start = round(half_max_indice_left);
    idx_end = round(ig);
    
    % зйщеб FWHM (аешк дфемс доевбм-фешййд)
    tau_TL = time_fs(idx_end) - time_fs(idx_start);
    t_start = time_fs(idx_start);
    t_end = time_fs(idx_end);
end

% =======================================================
% --- щмб 5: дцвъ ъецаеъ евшфйн ---
% =======================================================

disp('===================================================');
fprintf('Calculated Time-Limited Pulse Duration (FWHM): %.2f fs\n', tau_TL);
disp('===================================================');

figure;

subplot(2,1,1);
plot(wavelengths_nm, Y_spectrum, 'k', 'LineWidth', 2);
title('Extracted Spectral Intensity vs. Wavelength');
xlabel('Wavelength (nm)');
ylabel('Intensity (a.u.)');
grid on;

subplot(2,1,2);
plot(time_fs, I_time_norm, 'r', 'LineWidth', 2);
hold on;
plot([t_start t_end], [I_half I_half], 'b--', 'LineWidth', 1.5); % че FWHM
title('Time-Limited Pulse Intensity Profile');
xlabel('Time (fs)');
ylabel('Normalized Intensity');
legend('Pulse Intensity', ['FWHM = ' num2str(tau_TL, '%.2f') ' fs'], 'Location', 'best');
xlim([-5*tau_TL, 5*tau_TL]); % двбмъ дцйш миеез шмеерий
grid on;



% 
% 
% % =======================================================
% % --- щмб 2: дошд оошзб аешк вм (lambda) мошзб ъгш (omega) ---
% % =======================================================
% 
% % 1. зщб аъ дъгшйн (omega)
% omega_orig = 2 * pi * c ./ wavelengths_nm; 
% E_omega_magnitude_orig = sqrt(Y_spectrum);
% 
% % 2. оййп аъ дъгшйн (omega) еаъ дощштъ (E) бдъан
% [omega_orig_sorted, sort_idx] = sort(omega_orig);
% E_omega_magnitude_orig_sorted = E_omega_magnitude_orig(sort_idx);
% 
% % 3. бцт айришфемцйд тм дръерйн дооейрйн:
% % E_omega_TL_Interpolated = interp1(omega_orig_sorted, E_omega_magnitude_orig_sorted, omega_lin, 'pchip', 0);
% 
% % =======================================================
% % --- щмб 3: айришфемцйд егвйод азйгд (зйерй м-FFT) ---
% % =======================================================
% 
% % FFT гешщ щ-omega йдйд ошеез баефп азйг (мйрйашй).
% 
% % а. чбйтъ шщъ д-omega дзгщд (мйрйашйъ)
% N_points = 2^10; % осфш рчегеъ гвйод (тгйу зжчд щм 2 м-FFT йтйм)
% omega_min = min(omega_orig);
% omega_max = max(omega_orig);
% d_omega = (omega_max - omega_min) / (N_points - 1);
% 
% % цйш аеовд дзгщ едошеез мйрйашйъ
% omega_lin = linspace(omega_min, omega_max, N_points);
% 
% % б. айришфемцйд щм ощштъ дщгд мшщъ дмйрйашйъ   
% % (дщъощ б-interp1, афщш мдщъощ б-'spline' ае 'linear')
% E_omega_TL = interp1(omega_orig, E_omega_magnitude_orig, omega_lin, 'pchip', 0);
% % figure; plot(E_omega_TL)
% % TL - Time-Limited (фажд афс)
% 
% % =======================================================
% % --- щмб 4: ишрсфешн фешййд дфек (IFFT) езйщеб цйш джоп ---
% % =======================================================
% 
% % 1. зйщеб цйш джоп (Time Axis)
% % цтг джоп (dt) дфек миеез дъгшйн длемм
% T_range = 2 * pi / d_omega;
% dt = T_range / N_points;
% % цйш жоп ооешлж сбйб 0, бйзйгеъ фоиещрйеъ (fs)
% time_fs = (-N_points/2 : N_points/2 - 1) * dt;
% 
% % 2. бйцет IFFT
% % йщ мошлж аъ дръерйн (fftshift) мфрй д-IFFT.
% % аре ощъощйн б-ifft(ifftshift(E)) лгй мгоеъ IFT швйм оеошлж.
% E_time_uncentered = ifft(E_omega_TL);
% 
% % ошлеж дъецад
% E_time = fftshift(E_time_uncentered);
% 
% % 3. зйщеб тецоъ дфемс бжоп
% I_time = abs(E_time).^2; 
% 
% % =======================================================
% % --- щмб 5: оцйаъ ощк дфемс (FWHM) ---
% % =======================================================
% 
% % щйоещ бферчцйд оебрйъ 'fwhm' (ан жойрд бвшсд щмк) ае зйщеб йгрй:
% 
% % ршоем дтецод
% I_time_norm = I_time / max(I_time);
% 
% % оцйаъ айргчсйн отм зцй двебд (0.5)
% % half_max_indices = find(I_time_norm >= 0.5);
% [~,indx]= max(I_time_norm);
% half_max_indice_left = interp1(I_time_norm(1:indx),1:indx,0.5);
% half_max_indice_right = interp1(I_time_norm(indx:end),indx:length(I_time_norm),0.5);
% 
% if isempty(half_max_indice_left || half_max_indice_right)
%     tau_TL = NaN;
% else
%     % дайргчс дшащеп едазшеп
%     idx_start = floor(half_max_indice_left);
%     idx_end = ceil(half_max_indice_right);
% 
%     % зйщеб FWHM т"й зйсеш джорйн доъайойн (айришфемцйд йлемд мдйеъ огейчъ йеъш)
%     t_start = time_fs(idx_start);
%     t_end = time_fs(idx_end);
% 
%     tau_TL = t_end - t_start;
% end
% 
% % =======================================================
% % --- щмб 6: дцвъ ъецаеъ евшфйн ---
% % =======================================================
% 
% disp('===================================================');
% disp('Fourier-Limited Pulse Duration Calculation Results');
% disp('===================================================');
% fprintf('Calculated Time-Limited Pulse Duration (FWHM): %.2f fs\n', tau_TL);
% 
% figure;
% 
% subplot(2,1,1);
% plot(wavelengths_nm, Y_spectrum, 'k', 'LineWidth', 2);
% title('Extracted Spectral Intensity vs. Wavelength');
% xlabel('Wavelength (nm)');
% ylabel('Intensity (a.u.)');
% 
% subplot(2,1,2);
% plot(time_fs, I_time_norm, 'r', 'LineWidth', 2);
% hold on;
% plot([t_start t_end], [0.5 0.5], 'b--', 'LineWidth', 1.5); % че FWHM
% 
end