
NumberOfRandomSymbols = 50;

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
plot(Symbols, 'r+');
grid on;
axis([-2 2 -2 2]);
xlabel('Real');
ylabel('Imaginary');
title('Symbol Constellation');
print('-dpng', '~/study/university/semester7/diccom/symbol_constellation.png');

[ak, bk] = encoder(SymbolBits, Symbols, SymbolBitValues, BitsPerSymbol);
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
[q_i, q_q] = matched_demodulate( ModulatedSymbolBits , 1, A_c, omega_c,  0, Ts, f_s);

% Decide and decode
DecodedBits = decoder(MLLDecision([q_i, q_q]', Symbols),Symbols,SymbolBitVector);

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
stem(q_i);
grid on;
set(gca, 'XTick', [1:length(q_i)]);
set(gca, 'YTick', [-1 -1/sqrt(2) 0 1/sqrt(2) 1]);
set(gca, 'YTickLabel', {'-1','-1/sqrt(2)','0','1/sqrt(2)','1'});
xlabel('k (Symbol number)');
ylabel('q_k^i');
axis([0 length(q_i)+1 -1 1]);
title('Sampled in-phase component of match filter output');
subplot(2,1,2);
stem(q_q);
grid on;
set(gca, 'XTick', [1:length(q_q)]);
set(gca, 'YTick', [-1 -1/sqrt(2) 0 1/sqrt(2) 1]);
set(gca, 'YTickLabel', {'-1','-1/sqrt(2)','0','1/sqrt(2)','1'});
xlabel('k (Symbol number)');
ylabel('q_k^q');
axis([0 length(q_i)+1 -1 1]);
title('Sampled quadrature component of match filter output');
print('-dpng', '~/study/university/semester7/diccom/sampled_matched_filter_output.png');
