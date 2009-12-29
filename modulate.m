function [ s_M ] = modulate( ak, bk, Ts, omega_c, phi_0 , dt)
%modulate Modulates using (cos, -sin) basis the encoded input vector (Ak = ak + j bk)
%   Modulates the symbols (ak + jbk), with a time resolution of dt

NumberOfSymbols = length(ak);

s_M = []
for i = 1:NumberOfSymbols
    a = ak(i);
    b = bk(i);
    SymbolTimes = i*Ts:dt:(i+1)*Ts;
    samples = a*cos(omega_c * SymbolTimes + phi_0) - b*sin(omega_c * SymbolTimes + phi_0);
    s_M = [s_M, samples];
end

end

