
NumberOfRandomSymbols = 5000;

% Modulation parameters
Rb = 800;
M = 4;
Bch = 800;

f_c = 20*10^3;% carrier frequency
omega_c = 2*pi*f_c; 

f_s = 8*f_c; % sampling frequency = 4 * carrier freq. so that the fft will give the peak in the middle
P_c = 10; % transmission power
A_c = sqrt(2*P_c); % transmission amplitude


BitsPerSymbol = log2(M);
Rs = Rb / BitsPerSymbol;
Ts = 1 / Rs;

SqrtHalf = 1/sqrt(2);

Symbols = [SqrtHalf + 1j * SqrtHalf, 
          -SqrtHalf + 1j * SqrtHalf, 
          -SqrtHalf - 1j * SqrtHalf, 
           SqrtHalf - 1j * SqrtHalf];
% Bits as values
SymbolBitValues = [0,1,3,2];
% Bits corresponding to symbols in a vector (used for decoder map)
SymbolBitVector = [0, 0, 0, 1, 1, 1, 1, 0];

%--------------------------------------------------


% Small data set
SymbolBits = [0, 0, 1, 1, 1, 0, 0, 1];

% Random bits data set
NumberOfRandomBits = NumberOfRandomSymbols*log2(M);
RandomBits = arrayfun(@(x) (x > 0.5), rand(NumberOfRandomBits, 1));

NumberOfOnes = sum(RandomBits);

% Print so we can verify that the distribution is uniform
disp('Ratio of ones / zeros = ');
disp(NumberOfOnes / NumberOfRandomBits);

% Plot the constellation
subplot(1,1,1);
plot(Symbols, 'r+', 'MarkerSize',17);
grid on;
axis([-1.5 1.5 -1.5 1.5]);
rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'LineWidth',1,'LineStyle',':');
daspect([1,1,1]);
xlabel('Real');
ylabel('Imaginary');
title('Symbol Constellation');
print('-dpng', '~/study/university/semester7/diccom/symbol_constellation.png');

[ak, bk] = encoder(SymbolBits, Symbols, SymbolBitValues, BitsPerSymbol);
SmallDataSetSymbols = ak + 1j * bk;
[ModulatedSymbolBits, t] = modulate(ak, bk, Ts, omega_c, 0, A_c, f_s);
subplot(1,1,1);
plot(t, ModulatedSymbolBits);
grid on;
xlabel('time');
ylabel('s_M(t)');
title('Modulated small data set in time domain');
print('-dpng', '~/study/university/semester7/diccom/modulated_small_dataset.png');


% Modulate the RandomBits data set
[akr, bkr] = encoder(RandomBits, Symbols, SymbolBitValues, BitsPerSymbol);
[ModulatedRandomBits, t]= modulate(akr, bkr, Ts, omega_c, 0, A_c, f_s);
N = length(ModulatedRandomBits);
FFTResult = fftshift(fft(ModulatedRandomBits) / N);
FreqResp = abs(FFTResult);
PhaseResp = angle(FFTResult);
% Fix the fft result so the we have negative frequencies on the left
frequencies = f_s / N * [-N/2 : N/2-1];

subplot(2,1,1);
plot(frequencies, FreqResp);
grid on;
xlabel('frequency');
ylabel('|S_M(f)|');
title('Magnitude of Modulated random data set in frequency domain');

subplot(2,1,2);
plot(frequencies, PhaseResp);
grid on;
xlabel('frequency');
ylabel('|S_M(f)|');
title('Phase of Modulated random data set in frequency domain');
print('-dpng', '~/study/university/semester7/diccom/modulated_random_dataset_fft.png');

%--------------------------------------------------------------------------
% Demodulate and match
[qk_i, qk_q] = matched_demodulate( ModulatedSymbolBits , 1, A_c, omega_c,  0, Ts, f_s);

% Decide and decode
DecodedBits = decoder(MLLDecision([qk_i, qk_q]', Symbols),Symbols,SymbolBitVector);

display('The decoded bits are:');
DecodedBits

% Show the transmit filter output
[s_di, s_dq ,t] = transmit_filter(ak, bk, Ts, f_s);
subplot(2,1,1);
plot(t, s_di);
grid on;
xlabel('time');
ylabel('s_d_i(t)');
title('In-phase component of transmit filter output');
subplot(2,1,2);
grid on;
plot(t, s_dq);
ylabel('s_d_q(t)');
title('Quadrature component of transmit filter output');
print('-dpng', '~/study/university/semester7/diccom/transmit_filter_output.png');


% Show the output of the matched filter
subplot(2,1,1);
stem(qk_i);
grid on;
set(gca, 'XTick', [1:length(qk_i)]);
set(gca, 'YTick', [-1 -1/sqrt(2) 0 1/sqrt(2) 1]);
set(gca, 'YTickLabel', {'-1','-1/sqrt(2)','0','1/sqrt(2)','1'});
xlabel('k (Symbol number)');
ylabel('q_k^i');
axis([0 length(qk_i)+1 -1 1]);
title('Sampled in-phase component of match filter output');
subplot(2,1,2);
stem(qk_q);
grid on;
set(gca, 'XTick', [1:length(qk_q)]);
set(gca, 'YTick', [-1 -1/sqrt(2) 0 1/sqrt(2) 1]);
set(gca, 'YTickLabel', {'-1','-1/sqrt(2)','0','1/sqrt(2)','1'});
xlabel('k (Symbol number)');
ylabel('q_k^q');
axis([0 length(qk_q)+1 -1 1]);
title('Sampled quadrature component of match filter output');
print('-dpng', '~/study/university/semester7/diccom/sampled_matched_filter_output.png');


% Show the matched filter outputs before sampling
[q_i, q_q, t] = matched_demodulate_nosampling( ModulatedSymbolBits , 1, A_c, omega_c,  0, Ts, f_s);
subplot(2,1,1);
plot(t, q_i);
grid on;
xlabel('time');
ylabel('q_i(t)');
title('In-phase component of matched filter output');
subplot(2,1,2);
plot(t, q_q);
grid on;
ylabel('q_q(t)');
title('Quadrature component of matched filter output');
print('-dpng', '~/study/university/semester7/diccom/matched_filter_output.png');

% 2 - no phase difference
% Constellation of received symbols before decision device
ReceivedSymbols = qk_i + 1j * qk_q;

display('Received symbols:');
display(ReceivedSymbols);

display('Distance of received symbols from transmitted symbols:');
display(arrayfun(@abs, SmallDataSetSymbols-ReceivedSymbols));

subplot(1,1,1);
plot(qk_i, qk_q, 'b*', 'MarkerSize', 14);
rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'LineWidth',1,'LineStyle',':');
daspect([1,1,1]);
grid on;
axis([-1.5 1.5 -1.5 1.5]);
xlabel('Real');
ylabel('Imaginary');
title('Received Symbols Constellation (ideal synchronized phase)');
print('-dpng', '~/study/university/semester7/diccom/recv_symbol_constellation.png');


% Compare constellations
subplot(1,1,1);
plot(real(Symbols), imag(Symbols), 'r+', qk_i, qk_q, 'b*', 'MarkerSize', 12);
rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'LineWidth',1,'LineStyle',':');
daspect([1,1,1]);
grid on;
axis([0.65 0.75 0.65 0.75]);
xlabel('Real');
ylabel('Imaginary');
legend('Transmitted','Received');
title('Transmitted and Received Symbols Constellation (ideal synchronized phase)');
print('-dpng', '~/study/university/semester7/diccom/trans_recv_symbol_constellation.png');



% 3 - phase difference

% Demodulate and match
[qk3_i, qk3_q] = matched_demodulate( ModulatedSymbolBits , 1, A_c, omega_c,  pi/6, Ts, f_s);

% Decide and decode
DecodedBits3 = decoder(MLLDecision([qk3_i, qk3_q]', Symbols),Symbols,SymbolBitVector);


% Constellation of received symbols before decision device
ReceivedSymbols3 = qk3_i + 1j * qk3_q;

display('Received symbols with phase difference pi/6:');
display(ReceivedSymbols3);

display('Distance of received symbols from transmitted symbols, with phase difference of pi/6:');
display(arrayfun(@abs, SmallDataSetSymbols - ReceivedSymbols3));

subplot(1,1,1);
plot(qk3_i, qk3_q, 'b*', 'MarkerSize', 14);
rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'LineWidth',1,'LineStyle',':');
daspect([1,1,1]);
grid on;
axis([-1.5 1.5 -1.5 1.5]);
xlabel('Real');
ylabel('Imaginary');
title('Received Symbols Constellation, for phase difference of pi/6');
print('-dpng', '~/study/university/semester7/diccom/recv_symbol_constellation_phase.png');


% Compare constellations
subplot(1,1,1);
plot(real(Symbols), imag(Symbols), 'r+', qk3_i, qk3_q, 'b*', 'MarkerSize', 12);
rectangle('Position',[-1,-1,2,2],'Curvature',[1,1],'LineWidth',1,'LineStyle',':');
daspect([1,1,1]);
grid on;
axis([-1.5 1.5 -1.5 1.5]);
xlabel('Real');
ylabel('Imaginary');
legend('Transmitted','Received');
title('Transmitted and Received Symbols Constellation with phase difference pi/6');
print('-dpng', '~/study/university/semester7/diccom/trans_recv_symbol_constellation_phase.png');


% Channel with noise - synced phase
%----------------------------------
% We want the noise centered at fc = 20khz (40pi*10^3 rad/sec)
% and with a bandwith of Bch = 800hz, so f1 = 19.6khz, f2 = 20.4khz.
% The sampling frequency is fs = 8fc.
% The filter accepts the bandwith values as factors of fs/2 (for a value of 1 we get fs/2, 0.5 -> fs/4,
% etc.)
% In our case fs/2 = 4fc = 1 [normalized], so f / 4fc = [normalized value]

f1 = 19.6 / 80;
f2 = 20.4 / 80;
ChannelWindowFreqs = [f1,f2];
ChannelWindowFilter = fir1(70,ChannelWindowFreqs);

% how we picked the number 70:
% arrayfun(@(x) (std(filter(fir1(x,ChannelWindowFreqs), 1, noise))^2/800), [20:10:70])
% Gives a list of N_0 values after filtering. We picked the one that was closest to the intended one (the
% one before filtering).

gamma_d_max = 10*log10( 2*erfcinv(1e-3)^2 );
gamma_d_min = 10*log10( 2*erfcinv(0.2)^2 );

gamma_d_vec = gamma_d_min : 2 : gamma_d_max;
gamma_d_linear_vec = 10.^(gamma_d_vec./10);
snr_bit_linear_vec = gamma_d_linear_vec/2;
snr_bit_db_vec = 10*log10(snr_bit_linear_vec);

P_r = P_c;
NoiseVariances = P_r*4*f_c ./ (Rs*gamma_d_linear_vec);


BERsAverage = [];
BERsStdDev = [];
% For each gamma_d
for i = 1:length(gamma_d_vec)
  BERs = [];
  % Get average and std deviation of several attempts
  for j = 1:5
    % The number of noise samples should equal the number of modulated data samples = length(ModulatedRandomBits)
    noise = sqrt(NoiseVariances(i)) .* randn(1, length(ModulatedRandomBits));
    NoiseSamples = filter(ChannelWindowFilter, 1, noise);
    ReceivedSignal = ModulatedRandomBits + NoiseSamples;

    % Demodulate and match
    [qkn_i, qkn_q] = matched_demodulate( ReceivedSignal , 1, A_c, omega_c,  0, Ts, f_s);
    % Decide and decode
    DecodedNoisyBits = decoder(MLLDecision([qkn_i, qkn_q]', Symbols),Symbols,SymbolBitVector);

    % Count errors
    BERs(j) = sum(abs(DecodedNoisyBits' - RandomBits))./(length(DecodedNoisyBits));
  end
  BERsAverage(i) =  mean(BERs);
  BERsStdDev(i) = std(BERs);
end

AnalyticBER = erfc(sqrt(gamma_d_linear_vec/2))/2;

subplot(1,1,1);
%loglog(gamma_d_linear_vec/2, BitErrorRates, gamma_d_linear_vec/2, );
hold off;
semilogy(snr_bit_db_vec, AnalyticBER, 'b-');
hold;
errorbar(snr_bit_db_vec, BERsAverage, BERsStdDev, 'r');

legend('Analytic', 'Simulation');
grid on;
xlabel('SND_d^{(bit)} [dB]');
ylabel('P_{err}^{(bit)}');
title('Bit Error Rate for noisy channel with synchronized phase');
%print('-dpng', '~/study/university/semester7/diccom/BER_noise.png');
hold off;



% % Channel with noise - unsynced phase
% %----------------------------------
subplot(1,1,1);
hold off;
semilogy(snr_bit_db_vec, AnalyticBER, 'b-');
hold;

PhaseBERsAverage = [];
PhaseBERsStdDev = [];
Phases = [5 10];
Colors = 'rgy';
% Perform the same calculation, once for each phase difference
for k = 1:length(Phases)
    PhaseDifference = 2*pi/360 .* Phases(k);
    % For each gamma_d
    for i = 1:length(gamma_d_vec)
        BERs = [];
        % Get average and std deviation of several attemps
        for j = 1:5
            % The number of noise samples should equal the number of modulated data samples = length(ModulatedRandomBits)
            noise = sqrt(NoiseVariances(i)) .* randn(1, length(ModulatedRandomBits));
            NoiseSamples = filter(ChannelWindowFilter, 1, noise);
            ReceivedSignal = ModulatedRandomBits + NoiseSamples;

            % Demodulate and match
            [qkn_i, qkn_q] = matched_demodulate( ReceivedSignal , 1, A_c, omega_c,  PhaseDifference, Ts, f_s);
            % Decide and decode
            DecodedNoisyBits = decoder(MLLDecision([qkn_i, qkn_q]', Symbols),Symbols,SymbolBitVector);

            % Count errors
            BERs(j) = sum(abs(DecodedNoisyBits' - RandomBits))/(length(DecodedNoisyBits));
        end
        PhaseBERsAverage(k,i) =  mean(BERs);
        PhaseBERsStdDev(k,i) = std(BERs);
    end
    errorbar(snr_bit_db_vec, PhaseBERsAverage(k,:), PhaseBERsStdDev(k,:), Colors(k));
end

legend('Analytic', 'Simulation - phase 5^o', 'Simulation - phase 10^o');
grid on;
xlabel('SND_d^{(bit)} [dB]');
ylabel('P_{err}^{(bit)}');
title('Bit Error Rate for noisy channel with unsynchronized phase');
%print('-dpng', '~/study/university/semester7/diccom/BER_noise_phase.png');
hold off;


% ----------------------------------------------------------------
