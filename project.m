
NumberOfRandomSymbols = 5;

% Modulation parameters
Rb = 800;
M = 4;
Bch = 800;


BitsPerSymbol = log2(M);
Rs = Rb / BitsPerSymbol;

SqrtHalf = 1/sqrt(2);

Symbols = [SqrtHalf + 1j * SqrtHalf, 
          -SqrtHalf + 1j * SqrtHalf, 
          -SqrtHalf - 1j * SqrtHalf, 
           SqrtHalf - 1j * SqrtHalf];
SymbolBits = [0,1,3,2];

% Random bits
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

[ak, bk] = encoder(RandomBits, Symbols, SymbolBits, BitsPerSymbol);

