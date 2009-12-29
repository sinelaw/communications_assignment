function [ q_ki, q_kq ] = matched_demodulate( s_r, A_c, omega_c,  psi, T_s, SamplesPerSecond)
%matched_demodulate Demodulates the 4QAM signal in s_r, using carrier freq.
%omega_c, LO amplitude A_c and phase offset psi.
%   Returns vector of pairs (2 x N) q_k = [[q1_ki; q1_kq], [q2_ki; q2_kq], ...]

dt = 1 / SamplesPerSecond;
NumberOfSamples = length(s_r);

SamplesPerSymbol = SamplesPerSecond * T_s;
NumberOfSymbols = NumberOfSamples / SamplesPerSymbol;

t = 0:dt:(NumberOfSamples*dt); % calculate time sample values
t(:,end) = []; % Erase last time value (belongs to next symbol)

% Demodulate
x_di =  A_c * cos(omega_c*t+psi) .* s_r;
plot(x_di);
x_dq = -A_c * sin(omega_c*t+psi) .* s_r;

% Receive filter (matched filter)
q_ki = zeros(NumberOfSymbols,1);
q_kq = zeros(NumberOfSymbols,1);
for i = 1:NumberOfSymbols
    % Sum intervals of length Ts, skipping to the interval of symbol number
    % i at each step
    k = i-1;
    lower_bound = k*SamplesPerSymbol+1;
    upper_bound = (k+1)*SamplesPerSymbol;
    TimeInterval = t(lower_bound : upper_bound);
    q_ki(i) = trapz( TimeInterval, x_di( lower_bound : upper_bound ) );
    q_kq(i) = trapz( TimeInterval, x_dq( lower_bound : upper_bound ) );
end

end

