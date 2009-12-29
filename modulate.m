function [ s_M, t ] = modulate( ak, bk, Ts, omega_c, phi_0 , A_c, SamplesPerSecond)
%modulate Modulates using (cos, -sin) basis the encoded input vector (Ak = ak + j bk)
%   Modulates the symbols (ak + jbk), with a time resolution of dt
%     omega_c = carrier frequency
%     phi_0   = carrier constant phase
%     A_c     = carrier amplitude

NumberOfSymbols = length(ak);
dt = 1 / SamplesPerSecond

% Pre-allocate arrays for better performance
s_M = zeros(1,Ts*(NumberOfSymbols)*SamplesPerSecond);
t = 0:dt:(Ts*(NumberOfSymbols)); % calculate time sample values
t(:,end) = []; % Erase last time value (belongs to next symbol)

% Modulation
for i = 1:NumberOfSymbols
    a = ak(i);
    b = bk(i);
    SymbolTimes = t(1,(i-1)*Ts*SamplesPerSecond+1: i*Ts*SamplesPerSecond);
    samples = A_c * ( a*cos(omega_c * SymbolTimes + phi_0) - b*sin(omega_c * SymbolTimes + phi_0));
    % Write the samples into big array
    s_M(1,(i-1)*Ts*SamplesPerSecond+1: i*Ts*SamplesPerSecond) = samples;
end

end

