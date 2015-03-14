% A simple RLS-LMS compression test.
% Uses sines embedded in noise.
% By Nicholas Esterer

clear;
N = 10000;  % length of signal
m = 32;     % number of sine tones to synthesize
A = 1;      % amplitude of sine tones
w0 = pi/64; % fundamental frequency
SNR = 18;   % signal to noise ratio in dB (base 10)
rho = sqrt((m*A^2)/(2*10^(SNR/10))); % standard deviation of noise 
lambda = 0.9999; % forgetting factor for RLS
p_rls = 4;       % RLS filter order
p_lms = 128;     % LMS filter order
be = 0.05;       % NLMS gradient coefficient

phi = rand(m,1)*2*pi; % phases of sinusoids
n = [1:N] - 1;        % sample indices
v = randn(1,N)*rho^2; % noise
x = sum(A*cos((1:m)'*n*w0+phi*ones(1,N)),1); % sum of sines
d = x + v; % sines plus noise

figure(1);

[A_,U] = rlspyk_analysis(d,p_rls,lambda);
subplot(2,2,1);
plot(n,A_(2:end,:));
title('RLS Prediction Coefficients (Analysis)');
xlabel('Sample number / Iteration');

[A_,E] = nlms_lp_analysis(U,be,p_lms);
subplot(2,2,3);
plot(n,A_);
title('NLMS Prediction Coefficients (Analysis)');
xlabel('Sample number / Iteration');

[A_,U] = nlms_lp_synthesis(E,be,p_lms);
subplot(2,2,4);
plot(n,A_);
title('NLMS Prediction Coefficients (Synthesis)');
xlabel('Sample number / Iteration');

[A_,D] = rlspyk_synthesis(U,p_rls,lambda);
subplot(2,2,2);
plot(n,A_(2:end,:));
title('RLS Prediction Coefficients (Synthesis)');
xlabel('Sample number / Iteration');

figure(2);

E_d_ratios = sum(reshape(E,100,100).^2)./sum(reshape(d,100,100).^2);
E_d_log_power = 6*log(E_d_ratios)/log(2);
plot(n,reshape(ones(100,1)*E_d_log_power,N,1));
title('Residual power (averaged over 100 sample windows)');
ylabel('dB Power');
xlabel('Sample number');
figure(3);

plot(n,6*log(D-d.')/log(2));
title('Reconstruction error');
ylabel('dB');
xlabel('Sample number');