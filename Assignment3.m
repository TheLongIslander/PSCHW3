% Homework 2 - Baseband Modulation

%% Step 1: Generate Random Bit Sequence
rng(208001821); % Replace RUID with your actual Rutgers ID
bb = randi([0 1], 1, 1000); % Generate 1000 random bits

%% Step 2: Generate Baseband Signals
T = 2; Ts = 0.5; fs = 1/Ts;
t = -T:Ts:T; % Time vector

p_t = sinc(t/T); % Sinc pulse
pc_t = (1 - cos(pi * t / T)) .* (t >= 0 & t <= T); % Given pulse

s_t = []; s_c_t = [];
for i = 1:10
    if bb(i) == 1
        s_t = [s_t, p_t]; s_c_t = [s_c_t, pc_t];
    else
        s_t = [s_t, -p_t]; s_c_t = [s_c_t, -pc_t];
    end
end
t_plot = 0:Ts:(length(s_t)-1)*Ts;

%% Step 3: Plot Signals
figure;
subplot(2,1,1); plot(t_plot, s_t, 'b');
xlabel('Time (s)'); ylabel('Amplitude'); title('Baseband Signal s(t)');
subplot(2,1,2); plot(t_plot, s_c_t, 'r');
xlabel('Time (s)'); ylabel('Amplitude'); title('Baseband Signal s_c(t)');

%% Step 4: Compute FFT
N = length(s_t);
S_f = fftshift(abs(fft(s_t, N))); S_c_f = fftshift(abs(fft(s_c_t, N)));
freq = linspace(-fs/2, fs/2, N);

figure;
subplot(2,1,1); semilogy(freq, S_f, 'b');
xlabel('Frequency (Hz)'); ylabel('|S(f)|'); title('FFT of s(t)');
subplot(2,1,2); semilogy(freq, S_c_f, 'r');
xlabel('Frequency (Hz)'); ylabel('|S_c(f)|'); title('FFT of s_c(t)');

%% Step 5: Compute Bandwidth and Energy
P_s = S_f.^2; P_sc = S_c_f.^2;
P_s_half = max(P_s) / 2; P_sc_half = max(P_sc) / 2;

% 3dB Bandwidth Calculation with Edge-Case Handling
if isempty(find(P_s >= P_s_half, 1, 'first')) || isempty(find(P_s >= P_s_half, 1, 'last'))
    bw_3dB_s = NaN; % Handle edge case
else
    bw_3dB_s = freq(find(P_s >= P_s_half, 1, 'last')) - freq(find(P_s >= P_s_half, 1, 'first'));
end

if isempty(find(P_sc >= P_sc_half, 1, 'first')) || isempty(find(P_sc >= P_sc_half, 1, 'last'))
    bw_3dB_sc = NaN;
else
    bw_3dB_sc = freq(find(P_sc >= P_sc_half, 1, 'last')) - freq(find(P_sc >= P_sc_half, 1, 'first'));
end

% 99% Bandwidth Calculation
P_s_total = sum(P_s); P_sc_total = sum(P_sc);
P_s_cumsum = cumsum(P_s) / P_s_total;
P_sc_cumsum = cumsum(P_sc) / P_sc_total;

bw_99_s = freq(find(P_s_cumsum >= 0.99, 1, 'first')) - freq(find(P_s_cumsum >= 0.01, 1, 'first'));
bw_99_sc = freq(find(P_sc_cumsum >= 0.99, 1, 'first')) - freq(find(P_sc_cumsum >= 0.01, 1, 'first'));

% Energy Calculation (Using Trapezoidal Integration)
E_s_f = trapz(freq, P_s);
E_sc_f = trapz(freq, P_sc);

% Display Results
disp(['3dB BW s(t): ', num2str(bw_3dB_s)]);
disp(['3dB BW s_c(t): ', num2str(bw_3dB_sc)]);
disp(['99% BW s(t): ', num2str(bw_99_s)]);
disp(['99% BW s_c(t): ', num2str(bw_99_sc)]);
disp(['Energy from FFT s(t): ', num2str(E_s_f)]);
disp(['Energy from FFT s_c(t): ', num2str(E_sc_f)]);
