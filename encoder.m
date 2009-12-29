function [ encoded ] = encoder( data , Symbols, SymbolBits , BitsPerSymbol)
%encoder Encode a stream of data
%   Encodes the 'data' row vector of bits, which is of the form [1,0,1,...]
%   The number of bits in 'data' must be even.
%   Symbols is a row vector of symbols, e.g. [1+1j, -1-1j, 0, 1]
%   SymbolBits is a row vector of encoding bit patterns (as numbers), e.g. [2,3,1,0]
%       each number will be converted to the symbol with the same index in
%       the Symbols array, i.e. 2 -> 1+1j, 3 -> -1-1j, etc.
%   Returns a N x BitsPerSymbol matrix, where N = number of encoded symbols


% convert from complex to [real, imag]
SymbolsCoeffs = [real(Symbols) , imag(Symbols)];
% Store in a cell array for Map
SymbolsCoeffsCel = mat2cell(SymbolsCoeffs, ones(size(SymbolsCoeffs,1),1)', 2);

% Create the map: binary code -> symbol
EncoderMap = containers.Map(num2cell(SymbolBits), SymbolsCoeffsCel);

% Encode the data
encoded = [];
for i = 1:BitsPerSymbol:length(data)
    bits = data(i: i + BitsPerSymbol - 1);
    % Convert from binary vector to a number
    value = polyval(bits, 2);
    % Lookup the value to find the symbol
    symbol = EncoderMap(value);
    
    encoded = [encoded ; symbol];
end

end

