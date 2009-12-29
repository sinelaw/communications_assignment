
NumberOfRandomSymbols = 5000;

% Modulation parameters
Rb = 800;
M = 4;
Bch = 800;
omega_c = 10000; % Depends on channel's parameters

BitsPerSymbol = log2(M);
Rs = Rb / BitsPerSymbol;
Ts = 1 / Rs;

SqrtHalf = 1/sqrt(2);

Symbols = [SqrtHalf + 1j * SqrtHalf, 
          -SqrtHalf + 1j * SqrtHalf, 
          -SqrtHalf - 1j * SqrtHalf, 
           SqrtHalf - 1j * SqrtHalf];
SymbolBitMap = [0,1,3,2];

% Small data set
SymbolBits = [0, 0, 0, 1, 1, 1, 1, 0];

% Random bits data set
NumberOfRandomBits = NumberOfRandomSymbols*log2(M);
RandomBits = arrayfun(@(x) (x > 0.5), rand(NumberOfRandomBits, 1));

NumberOfOnes = sum(RandomBits);

% Print so we can verify that the distribution is uniform
disp('Ratio of ones / zeros = ');
disp(NumberOfOnes / NumberOfRandomBits);

% Plot the constellation
plot(Symbols, 'r+');
grid on;
axis([-2 2 -2 2]);
xlabel('Real');
ylabel('Imaginary');
title('Symbol Constellation');
print('-dpng', '~/study/university/semester7/diccom/symbol_constellation.png');

[ak, bk] = encoder(SymbolBits, Symbols, SymbolBitMap, BitsPerSymbol);
[ModulatedSymbolBits, t] = modulate(ak, bk, Ts, omega_c, 0, 100);
plot(t, ModulatedSymbolBits);
xlabel('time');
ylabel('s_M(t)');
title('Modulated small data set in time domain');
print('-dpng', '~/study/university/semester7/diccom/modulated_small_dataset.png');


[ak, bk] = encoder(RandomBits, Symbols, SymbolBitMap, BitsPerSymbol);
[ModulatedRandomBits, t]= modulate(ak, bk, Ts, omega_c, 0, 100);
% TODO calculate f x-axis
plot(abs(fft(ModulatedRandomBits)));
xlabel('frequency');
ylabel('S_M(f)');
title('Modulated random data set in frequency domain');
print('-dpng', '~/study/university/semester7/diccom/modulated_random_dataset_fft.png');
