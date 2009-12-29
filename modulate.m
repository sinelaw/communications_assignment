function [ s_M, t ] = modulate( ak, bk, Ts, omega_c, phi_0 , A_c, SamplesPerSecond)
%modulate Modulates using (cos, -sin) basis the encoded input vector (Ak = ak + j bk)
%   Modulates the symbols (ak + jbk), with a time resolution of dt
%     omega_c = carrier frequency
%     phi_0   = carrier constant phase
%     A_c     = carrier amplitude

NumberOfSymbols = length(ak);
dt = 1 / SamplesPerSecond;

% Pre-allocate arrays for better performance
% using round to prevent matlab warning about non-integers (although the
% result must be an integer)
s_M = zeros(1,round(Ts*NumberOfSymbols*SamplesPerSecond));
t = 0:dt:Ts*NumberOfSymbols; % calculate time sample values
t(:,end) = []; % Erase last time value (belongs to next symbol)

% Modulation
for i = 1:NumberOfSymbols
    a = ak(i);
    b = bk(i);
    
    % The bounds are really always integers, but floating point errors 
    % cause matlab to print warnings
    lower_bound = round((i-1)*Ts*SamplesPerSecond+1);
    upper_bound = round(i*Ts*SamplesPerSecond);
    
    SymbolTimes = t(1,lower_bound:upper_bound);
    samples = A_c * ( a*cos(omega_c * SymbolTimes + phi_0) - b*sin(omega_c * SymbolTimes + phi_0));
    % Write the samples into big array
    s_M(1,lower_bound:upper_bound) = samples;
end

end

