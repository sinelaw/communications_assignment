function [ b_r ] = decoder( symbols_data , Symbols, SymbolBitVectpr)
%decoder Decode a stream of symbols
%   Decodes the 'symbols_data' row vector of symbols, which is of the form [1+1j,-1-1j,...]
%   Symbols is a row vector of symbols, e.g. [1+1j, -1-1j, 0, 1]
%   SymbolBitVectpr is a row vector of encoding bit patterns (as numbers), e.g. [0,0,1,1,...]
%       each symbol will be converted to the bit pattern with the same index in
%       the Symbols array, i.e. 2 <- 1+1j, 3 <- -1-1j, etc.
%   Returns a row vector of bits


% Store in a cell array for Map
SymbolBitVectprCoeffsCel = mat2cell(SymbolBitVectpr, 1, 2*ones(length(SymbolBitVectpr)/2,1));

% Create the map: binary code -> symbol
decoderMap = containers.Map(num2cell(1:length(Symbols)), SymbolBitVectprCoeffsCel);

% Encode the data
b_r = [];
for i = 1:length(symbols_data)
    symbol = symbols_data(i);
    % Lookup the value to find the symbol
    bits = decoderMap(symbol);
    b_r = [b_r, bits];
end

end

