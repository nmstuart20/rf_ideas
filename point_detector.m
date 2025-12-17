% RF Signal Detection via Triangle Wave Identification
% Detects a constant bandwidth signal by convolving FFT magnitude with 
% a square wave, then identifying the resulting triangle wave

%% Parameters
fs = 1e6;                    % Sampling frequency (Hz)
N = 2048;                    % FFT size
signal_bw = 50e3;            % Known signal bandwidth (Hz)
freq_range = [100e3, 500e3]; % Frequency search range (Hz)

% Detection thresholds
peak_sharpness_threshold = 0.7;  % For pointiness detection (normalized)
correlation_threshold = 0.8;      % For matched filter detection (normalized)

%% Generate Test Signal (for demonstration)
% Remove this section and use your actual I&Q data
t = (0:N-1)/fs;
signal_freq = 250e3;  % Actual signal frequency (unknown in practice)
% Simulate bandlimited signal
iq_signal = randn(1,N) + 1j*randn(1,N);  % Complex baseband noise
% Add signal component
f = linspace(-fs/2, fs/2, N);
signal_mask = abs(f - signal_freq) < signal_bw/2;
IQ_fft = fft(iq_signal);
IQ_fft = IQ_fft .* fftshift(signal_mask)' * 10;  % Boost signal region
iq_signal = ifft(IQ_fft);

%% Step 1: Compute FFT of I&Q data
IQ_fft = fft(iq_signal, N);
magnitude_spectrum = abs(fftshift(IQ_fft));
freq_axis = linspace(-fs/2, fs/2, N);

%% Step 2: Create square wave template with known bandwidth
freq_resolution = fs / N;
bw_bins = round(signal_bw / freq_resolution);
square_template = zeros(1, N);
center_idx = floor(N/2) + 1;
half_bw_bins = floor(bw_bins/2);
square_template(center_idx-half_bw_bins:center_idx+half_bw_bins) = 1;

%% Step 3: Convolve with square wave
convolution_result = conv(magnitude_spectrum, square_template, 'same');
convolution_result = convolution_result / max(convolution_result);  % Normalize

%% Step 4: Detection Method 1 - Peak Sharpness (Pointiness)
% Find peaks in convolution result
[pks, locs] = findpeaks(convolution_result, 'MinPeakProminence', 0.3);

% Calculate sharpness for each peak
sharpness_scores = zeros(size(pks));
for i = 1:length(locs)
    idx = locs(i);
    % Extract region around peak
    window_size = min(bw_bins * 2, 50);  % Adaptive window
    idx_start = max(1, idx - window_size);
    idx_end = min(N, idx + window_size);
    peak_region = convolution_result(idx_start:idx_end);
    
    % Calculate sharpness metric: ratio of peak value to width at half max
    half_max = pks(i) / 2;
    above_half = peak_region > half_max;
    width_at_half = sum(above_half);
    
    % Sharpness: higher value = sharper peak (triangle-like)
    if width_at_half > 0
        sharpness_scores(i) = pks(i) / (width_at_half * freq_resolution);
    end
end

% Normalize sharpness scores
if ~isempty(sharpness_scores)
    sharpness_scores = sharpness_scores / max(sharpness_scores);
end

% Detect signals based on sharpness
detected_method1 = sharpness_scores > peak_sharpness_threshold;
detected_freqs_method1 = freq_axis(locs(detected_method1));

%% Step 5: Detection Method 2 - Matched Filter (Triangle Template)
% Create ideal triangle wave template
triangle_template = zeros(1, 2*bw_bins+1);
for i = 1:length(triangle_template)
    if i <= bw_bins+1
        triangle_template(i) = i;  % Rising edge
    else
        triangle_template(i) = 2*bw_bins + 2 - i;  % Falling edge
    end
end
triangle_template = triangle_template / max(triangle_template);  % Normalize

% Correlate convolution result with triangle template
correlation = conv(convolution_result, triangle_template, 'same');
correlation = correlation / max(abs(correlation));  % Normalize

% Find peaks in correlation
[corr_pks, corr_locs] = findpeaks(abs(correlation), 'MinPeakProminence', 0.2);

% Detect signals based on correlation threshold
detected_method2 = corr_pks > correlation_threshold;
detected_freqs_method2 = freq_axis(corr_locs(detected_method2));

%% Step 6: Combine Detection Results
% A detection is confirmed if either method detects it
all_detected_freqs = unique([detected_freqs_method1, detected_freqs_method2]);

%% Display Results
fprintf('\n=== RF Signal Detection Results ===\n');
fprintf('Search Range: %.1f kHz to %.1f kHz\n', freq_range(1)/1e3, freq_range(2)/1e3);
fprintf('Known Signal Bandwidth: %.1f kHz\n', signal_bw/1e3);
fprintf('\nMethod 1 (Peak Sharpness): %d signal(s) detected\n', sum(detected_method1));
if any(detected_method1)
    for i = 1:length(detected_freqs_method1)
        fprintf('  - %.3f kHz (sharpness: %.3f)\n', ...
            detected_freqs_method1(i)/1e3, sharpness_scores(find(locs == locs(detected_method1), 1)));
    end
end

fprintf('\nMethod 2 (Matched Filter): %d signal(s) detected\n', sum(detected_method2));
if any(detected_method2)
    for i = 1:length(detected_freqs_method2)
        fprintf('  - %.3f kHz (correlation: %.3f)\n', ...
            detected_freqs_method2(i)/1e3, corr_pks(detected_method2));
    end
end

fprintf('\nCombined Detection: %d unique signal(s) detected\n', length(all_detected_freqs));
if ~isempty(all_detected_freqs)
    fprintf('Detected at frequencies:\n');
    for i = 1:length(all_detected_freqs)
        fprintf('  - %.3f kHz\n', all_detected_freqs(i)/1e3);
    end
end

%% Visualization
figure('Position', [100, 100, 1200, 800]);

% Plot 1: Original FFT Magnitude
subplot(4,1,1);
plot(freq_axis/1e3, magnitude_spectrum, 'b', 'LineWidth', 1);
xlabel('Frequency (kHz)');
ylabel('Magnitude');
title('FFT Magnitude Spectrum');
grid on;
xlim([freq_range(1)/1e3, freq_range(2)/1e3]);

% Plot 2: Convolution Result (should show triangle)
subplot(4,1,2);
plot(freq_axis/1e3, convolution_result, 'r', 'LineWidth', 1.5);
hold on;
if ~isempty(locs)
    plot(freq_axis(locs)/1e3, pks, 'ko', 'MarkerSize', 8, 'LineWidth', 2);
end
xlabel('Frequency (kHz)');
ylabel('Normalized Amplitude');
title('Convolution Result (Square * Spectrum) - Triangle Wave Expected');
grid on;
xlim([freq_range(1)/1e3, freq_range(2)/1e3]);

% Plot 3: Sharpness Analysis
subplot(4,1,3);
if ~isempty(locs)
    stem(freq_axis(locs)/1e3, sharpness_scores, 'filled');
    hold on;
    yline(peak_sharpness_threshold, 'r--', 'Threshold', 'LineWidth', 2);
    if any(detected_method1)
        plot(detected_freqs_method1/1e3, sharpness_scores(detected_method1), ...
            'go', 'MarkerSize', 12, 'LineWidth', 3);
    end
end
xlabel('Frequency (kHz)');
ylabel('Sharpness Score');
title('Method 1: Peak Sharpness Detection');
grid on;
xlim([freq_range(1)/1e3, freq_range(2)/1e3]);

% Plot 4: Matched Filter Correlation
subplot(4,1,4);
plot(freq_axis/1e3, abs(correlation), 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
hold on;
if ~isempty(corr_locs)
    plot(freq_axis(corr_locs)/1e3, corr_pks, 'ko', 'MarkerSize', 8, 'LineWidth', 2);
end
yline(correlation_threshold, 'r--', 'Threshold', 'LineWidth', 2);
if any(detected_method2)
    plot(detected_freqs_method2/1e3, corr_pks(detected_method2), ...
        'go', 'MarkerSize', 12, 'LineWidth', 3);
end
xlabel('Frequency (kHz)');
ylabel('Correlation');
title('Method 2: Matched Filter (Triangle Template)');
grid on;
xlim([freq_range(1)/1e3, freq_range(2)/1e3]);

%% Optional: Zoom in on detected signals
if ~isempty(all_detected_freqs)
    figure('Position', [150, 150, 1000, 600]);
    for i = 1:min(length(all_detected_freqs), 3)  % Show up to 3 detections
        subplot(min(length(all_detected_freqs), 3), 1, i);
        det_freq = all_detected_freqs(i);
        [~, center_idx] = min(abs(freq_axis - det_freq));
        zoom_range = round(bw_bins * 3);
        idx_range = max(1, center_idx-zoom_range):min(N, center_idx+zoom_range);
        
        plot(freq_axis(idx_range)/1e3, convolution_result(idx_range), 'r', 'LineWidth', 2);
        hold on;
        plot(det_freq/1e3, convolution_result(center_idx), 'go', ...
            'MarkerSize', 12, 'LineWidth', 3);
        xlabel('Frequency (kHz)');
        ylabel('Amplitude');
        title(sprintf('Detected Signal at %.3f kHz', det_freq/1e3));
        grid on;
    end
end
