function [ symbols, phase, phase_error, t ] = DecisionFeedbackLoopReceiver( s_r, A_r, A_0, omega_c, const_phase_offset, Ts, ...
                                                                 phase_variance, SNRbit, Symbols, SamplesPerSecond)
%DecisionFeedbackLoopReceiver Receives a modulated data stream in s_r and uses a DFL receiver to decode the data.
%   Returns vector of pairs (2 x N) q_k = [[q1_ki; q1_kq], [q2_ki; q2_kq], ...]

dt = 1 / SamplesPerSecond;
NumberOfSamples = length(s_r);

SamplesPerSymbol = SamplesPerSecond * Ts;
NumberOfSymbols = NumberOfSamples / SamplesPerSymbol;

t = 0:dt:(NumberOfSamples*dt); % calculate time sample values
t(:,end) = []; % Erase last time value (belongs to next symbol)

SingleTimeInterval = 0:dt:Ts;
SingleTimeInterval(:,end) = [];

% DFL Parameters
Rs = 1/Ts;
Pr = A_r^2 / 2;
N0 = Pr / (2*Rs*SNRbit);
Kpd = A_r * A_0 / 2;
tau1 = 0.5;
tau2 = 7.12e-3;
Beq = Kpd^2*phase_variance/N0;
K = 96.6; %4*Beq;
Kvco = K/Kpd;
% ---------------------

LoopFilter = tf([tau2,1],[tau1,1]); %, SamplesPerSecond);
VCOFilter = tf([Kvco],[1,0]); %, SamplesPerSecond);

q_ki = zeros(NumberOfSymbols,1);
q_kq = zeros(NumberOfSymbols,1);

EstimatedPhases = zeros(1,SamplesPerSymbol);
ak_i = -cos(pi/4);
ak_q = -sin(pi/4);

symbols = zeros(2, NumberOfSymbols);
phase_error = zeros(1, NumberOfSymbols*SamplesPerSymbol);
phase = zeros(1, NumberOfSymbols*SamplesPerSymbol);
for i = 1:NumberOfSymbols
    % i at each step
    k = i-1;
    lower_bound = k*SamplesPerSymbol+1;
    upper_bound = (k+1)*SamplesPerSymbol;
    TimeInterval = t(1, lower_bound : upper_bound);
    
    [s1a,s2a] = PhaseDetector(s_r(1,lower_bound:upper_bound), EstimatedPhases, A_0, omega_c, const_phase_offset, TimeInterval);
    
    s1b = filter(lowpass_1000, s1a) .* ak_q;
    s2b = filter(lowpass_1000, s2a) .* ak_i;

    % Phase error
    epsilon = filter(lowpass_1000, s2b - s1b);
    phase_error(1,lower_bound:upper_bound) = epsilon;

    LoopFilterOutput = lsim(LoopFilter, epsilon', SingleTimeInterval);
    NextEstimatedPhases = lsim(VCOFilter, LoopFilterOutput, SingleTimeInterval)';
    phase(1,lower_bound:upper_bound) = NextEstimatedPhases;

    % Demodulate the current symbol
    [qk_i, qk_q] = matched_demodulate( s_r(1,lower_bound:upper_bound) , A_0, A_r, omega_c, ...
                                       const_phase_offset, EstimatedPhases, Ts, SamplesPerSecond);
    [sym_index, symbol] = MLLDecision([qk_i, qk_q]', Symbols);
    
    ak_i = real(symbol);
    ak_q = imag(symbol);
    symbols(1, i) = qk_i;
    symbols(2, i) = qk_q;
    EstimatedPhases = NextEstimatedPhases;
end

end

