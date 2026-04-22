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

function zip_to_barcode(zipPath, outputImagePath, moduleSize)
% ZIP_TO_BARCODE Encode a zip file into a printable 2D barcode image.
%
%   zip_to_barcode(zipPath, outputImagePath, moduleSize)
%
%   Inputs:
%     zipPath         - path to the .zip file to encode
%     outputImagePath - path for the output PNG (e.g. ‘barcode.png’)
%     moduleSize      - (optional) pixels per data module (default 8)
%
%   The output is a 2D matrix-style barcode with:
%     - L-shaped finder pattern (solid left + bottom borders)
%     - Dotted timing pattern (top + right borders)
%     - A 64-bit header encoding the payload length in bytes
%     - CRC-32 checksum appended to the payload for integrity
%     - Quiet zone (white border) around the symbol
%
%   Uses ONLY built-in MATLAB functions (no toolboxes required).

```
if nargin < 3 || isempty(moduleSize)
    moduleSize = 8;
end

% ---- 1. Read the zip file as raw bytes -----------------------------
fid = fopen(zipPath, 'rb');
if fid < 0
    error('Could not open file: %s', zipPath);
end
payload = fread(fid, inf, 'uint8=>uint8');
fclose(fid);
payload = payload(:);               % column vector
nBytes  = numel(payload);
fprintf('Read %d bytes from %s\n', nBytes, zipPath);

% ---- 2. Append CRC-32 of the payload -------------------------------
crc = crc32(payload);               % uint32
crcBytes = typecast(swapbytes(crc), 'uint8');   % big-endian 4 bytes
payloadWithCrc = [payload; crcBytes(:)];

% ---- 3. Build header: 8 bytes = payload length (big-endian) --------
lenBytes = typecast(swapbytes(uint64(nBytes)), 'uint8');
fullStream = [lenBytes(:); payloadWithCrc];

% ---- 4. Convert bytes to a bit stream (MSB first) ------------------
bits = byteStreamToBits(fullStream);
nBits = numel(bits);

% ---- 5. Decide grid size: square, with reserved borders ------------
% Reserve 1 row/col on all four sides for finder + timing patterns.
% So usable data area = (N-2) x (N-2) bits.
N = ceil(sqrt(nBits)) + 2;
while (N-2)^2 < nBits
    N = N + 1;
end
fprintf('Encoding %d bits in a %dx%d module grid\n', nBits, N, N);

% ---- 6. Build the module grid (0 = white, 1 = black) ---------------
grid = zeros(N, N, 'uint8');

% Solid finder pattern: left column and bottom row
grid(:, 1)   = 1;
grid(end, :) = 1;

% Timing pattern: alternating on top row and right column
topRow = mod(0:N-1, 2);             % 0,1,0,1,...
grid(1, :) = topRow;
grid(:, end) = mod((0:N-1)', 2);

% Corners fixups so finder L is clean
grid(1, 1) = 1;                     % top-left start of L? keep black
grid(end, end) = 1;                 % bottom-right corner black
grid(1, end)   = 1;                 % top-right anchor black

% Fill data area row-by-row, left-to-right, top-to-bottom
dataRows = 2:N-1;
dataCols = 2:N-1;
bitIdx = 1;
for r = dataRows
    for c = dataCols
        if bitIdx <= nBits
            grid(r, c) = bits(bitIdx);
            bitIdx = bitIdx + 1;
        else
            % Pad remaining modules with a fixed pattern (checker)
            grid(r, c) = mod(r + c, 2);
        end
    end
end

% ---- 7. Upscale to moduleSize pixels per module --------------------
img = kron(grid, ones(moduleSize, moduleSize, 'uint8'));
% Convert to black/white image: 0 = black ink, 255 = white paper
img = uint8(255 * (1 - img));

% ---- 8. Add quiet zone (white border) ------------------------------
quiet = 4 * moduleSize;
[h, w] = size(img);
framed = 255 * ones(h + 2*quiet, w + 2*quiet, 'uint8');
framed(quiet+1:quiet+h, quiet+1:quiet+w) = img;

% ---- 9. Save as PNG ------------------------------------------------
imwrite(framed, outputImagePath);
fprintf('Saved barcode image (%dx%d px) to %s\n', ...
        size(framed,2), size(framed,1), outputImagePath);
fprintf('Grid: %dx%d modules, moduleSize=%d px\n', N, N, moduleSize);
```

end

% ======================================================================
function bits = byteStreamToBits(bytes)
% Convert a uint8 column vector into a column vector of bits (MSB first).
bytes = uint8(bytes(:));
n = numel(bytes);
bits = zeros(8*n, 1, ‘uint8’);
for k = 0:7
bits(k+1:8:end) = bitget(bytes, 8-k);
end
end

% ======================================================================
function c = crc32(data)
% CRC32 Compute the IEEE 802.3 CRC-32 of a uint8 vector.
data = uint8(data(:));
persistent table
if isempty(table)
table = zeros(256, 1, ‘uint32’);
poly  = uint32(hex2dec(‘EDB88320’));
for i = 0:255
r = uint32(i);
for j = 1:8
if bitand(r, uint32(1))
r = bitxor(bitshift(r, -1), poly);
else
r = bitshift(r, -1);
end
end
table(i+1) = r;
end
end
c = uint32(hex2dec(‘FFFFFFFF’));
for k = 1:numel(data)
idx = bitand(bitxor(c, uint32(data(k))), uint32(255)) + 1;
c = bitxor(bitshift(c, -8), table(idx));
end
c = bitxor(c, uint32(hex2dec(‘FFFFFFFF’)));
end

function barcode_to_zip(imagePath, outputZipPath)
% BARCODE_TO_ZIP Decode a 2D barcode image (made by zip_to_barcode) back
% into its original zip file.
%
%   barcode_to_zip(imagePath, outputZipPath)
%
%   Inputs:
%     imagePath     - path to the PNG produced by zip_to_barcode
%     outputZipPath - path for the recovered .zip
%
%   Uses ONLY built-in MATLAB functions.

```
% ---- 1. Load image and binarize -----------------------------------
img = imread(imagePath);
if ndims(img) == 3
    img = uint8(mean(img, 3));
end
% Black modules in the image are 0, white are 255 -> invert
% so that data-modules = 1, paper = 0.
bw = img < 128;

% ---- 2. Crop off the quiet zone (white border) --------------------
colHas = any(bw, 1);
rowHas = any(bw, 2);
c1 = find(colHas, 1, 'first');  c2 = find(colHas, 1, 'last');
r1 = find(rowHas, 1, 'first');  r2 = find(rowHas, 1, 'last');
if isempty(c1) || isempty(r1)
    error('No barcode content found in image.');
end
bw = bw(r1:r2, c1:c2);
[H, W] = size(bw);

% ---- 3. Determine module size from the solid bottom finder bar ----
% The bottom row of the symbol is one module tall. Scan UP from the
% bottom at a middle column until we exit black. The left column
% is also one module wide but a scan through the data area can be
% inflated by an adjacent black data module, so we use it only as
% a sanity check and take the minimum.
midCol = round(W/2);
botRun = 0;
for y = H:-1:1
    if bw(y, midCol)
        botRun = botRun + 1;
    else
        break;
    end
end
midRow = round(H/2);
leftRun = 0;
for x = 1:W
    if bw(midRow, x)
        leftRun = leftRun + 1;
    else
        break;
    end
end
if botRun >= 1
    moduleSize = min(leftRun, botRun);
else
    moduleSize = leftRun;
end
if moduleSize < 1
    error('Could not detect module size.');
end
fprintf('Detected moduleSize = %d pixels\n', moduleSize);

% ---- 4. Compute grid dimensions and sample each module ------------
N = round(H / moduleSize);
if round(W / moduleSize) ~= N
    % Allow slight rounding mismatch; use the average
    N = round((H + W) / (2 * moduleSize));
end
fprintf('Detected grid size = %dx%d modules\n', N, N);

grid = zeros(N, N, 'uint8');
for r = 1:N
    for c = 1:N
        % Sample center pixel of each module
        y = round((r - 0.5) * moduleSize);
        x = round((c - 0.5) * moduleSize);
        y = max(1, min(H, y));
        x = max(1, min(W, x));
        grid(r, c) = bw(y, x);
    end
end

% ---- 5. Extract data-area bits (inner (N-2)x(N-2) block) ----------
dataArea = grid(2:N-1, 2:N-1);
bits = reshape(dataArea', [], 1);     % row-major, matches encoder

% ---- 6. Convert bits back into bytes ------------------------------
nFullBytes = floor(numel(bits) / 8);
bits = bits(1:nFullBytes*8);
bytes = bitsToByteStream(bits);

% ---- 7. Parse header: first 8 bytes = payload length --------------
if numel(bytes) < 12
    error('Decoded data is too short to contain header + CRC.');
end
lenBytes = bytes(1:8);
nBytes = double(swapbytes(typecast(uint8(lenBytes), 'uint64')));
fprintf('Header says payload length = %d bytes\n', nBytes);

if numel(bytes) < 8 + nBytes + 4
    error(['Decoded data shorter than header claims ' ...
           '(have %d, need %d).'], numel(bytes), 8 + nBytes + 4);
end

payload  = bytes(9:8+nBytes);
crcBytes = bytes(8+nBytes+1 : 8+nBytes+4);
storedCrc = swapbytes(typecast(uint8(crcBytes), 'uint32'));

% ---- 8. Verify CRC-32 ---------------------------------------------
computedCrc = crc32(payload);
if storedCrc ~= computedCrc
    warning(['CRC mismatch! stored=0x%08X computed=0x%08X. ' ...
             'File may be corrupted.'], storedCrc, computedCrc);
else
    fprintf('CRC-32 OK (0x%08X)\n', storedCrc);
end

% ---- 9. Write the recovered zip -----------------------------------
fid = fopen(outputZipPath, 'wb');
if fid < 0
    error('Could not open output file: %s', outputZipPath);
end
fwrite(fid, payload, 'uint8');
fclose(fid);
fprintf('Recovered zip written to %s (%d bytes)\n', ...
        outputZipPath, numel(payload));
```

end

% ======================================================================
function bytes = bitsToByteStream(bits)
bits = uint8(bits(:));
n = floor(numel(bits) / 8);
bytes = zeros(n, 1, ‘uint8’);
for k = 0:7
bytes = bitor(bytes, bitshift(bits(k+1:8:8*n), 7-k));
end
end

% ======================================================================
function c = crc32(data)
data = uint8(data(:));
persistent table
if isempty(table)
table = zeros(256, 1, ‘uint32’);
poly  = uint32(hex2dec(‘EDB88320’));
for i = 0:255
r = uint32(i);
for j = 1:8
if bitand(r, uint32(1))
r = bitxor(bitshift(r, -1), poly);
else
r = bitshift(r, -1);
end
end
table(i+1) = r;
end
end
c = uint32(hex2dec(‘FFFFFFFF’));
for k = 1:numel(data)
idx = bitand(bitxor(c, uint32(data(k))), uint32(255)) + 1;
c = bitxor(bitshift(c, -8), table(idx));
end
c = bitxor(c, uint32(hex2dec(‘FFFFFFFF’)));
end









