function [ q_ki, q_kq ] = matched_demodulate( s_r, A_0, A_r, omega_c,  constant_phase,  varying_phase, Ts, SamplesPerSecond)
%matched_demodulate Demodulates the 4QAM signal in s_r, using carrier freq.
%omega_c, LO amplitude A_r and phase offset constant_phase.
%   Returns vector of pairs (2 x N) q_k = [[q1_ki; q1_kq], [q2_ki; q2_kq], ...]

dt = 1 / SamplesPerSecond;
NumberOfSamples = length(s_r);

SamplesPerSymbol = SamplesPerSecond * Ts;
NumberOfSymbols = NumberOfSamples / SamplesPerSymbol;

t = 0:dt:(NumberOfSamples*dt); % calculate time sample values
t(:,end) = []; % Erase last time value (belongs to next symbol)

% Demodulate. We multiply by 2Rs/A_r = 2/(A_r *Ts) already here, it's for restoring the
% amplitudes to the transmitted value (of course disregarding channel
% deamplification)
x_di =  A_0 * cos(omega_c*t + constant_phase .+ varying_phase) .* s_r * 2/ (A_r*Ts);
x_dq = -A_0 * sin(omega_c*t + constant_phase .+ varying_phase) .* s_r * 2/ (A_r*Ts);

% Receive filter (matched filter)
q_ki = zeros(NumberOfSymbols,1);
q_kq = zeros(NumberOfSymbols,1);
for i = 1:NumberOfSymbols
    % Sum intervals of length Ts, skipping to the interval of symbol number
    % i at each step
    k = i-1;
    lower_bound = k*SamplesPerSymbol+1;
    upper_bound = (k+1)*SamplesPerSymbol;
    TimeInterval = t(1,lower_bound : upper_bound);
    
    q_ki(i) = trapz( TimeInterval, x_di( lower_bound : upper_bound ) );
    q_kq(i) = trapz( TimeInterval, x_dq( lower_bound : upper_bound ) );
end

end

