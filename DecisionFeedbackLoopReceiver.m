function [ q_ki, q_kq ] = DecisionFeedbackLoopReceiver( s_r, A_r, A_0, omega_c, const_phase_offset, Ts, ...
                                                  phase_variance, SNRbit, SamplesPerSecond)
%DecisionFeedbackLoopReceiver Receives a modulated data stream in s_r and uses a DFL receiver to decode the data.
%   Returns vector of pairs (2 x N) q_k = [[q1_ki; q1_kq], [q2_ki; q2_kq], ...]

dt = 1 / SamplesPerSecond;
NumberOfSamples = length(s_r);

SamplesPerSymbol = SamplesPerSecond * Ts;
NumberOfSymbols = NumberOfSamples / SamplesPerSymbol;

t = 0:dt:(NumberOfSamples*dt); % calculate time sample values
t(:,end) = []; % Erase last time value (belongs to next symbol)

% DFL Parameters
Rs = 1/Ts;
Pr = A_r^2 / 2;
N0 = Pr / (2*Rs*SNRbit);
Kpd = A_r * A_0 / 2;
tau1 = 0.5;
tau2 = 0;
Beq = Kpd^2*phase_variance/N0;
K = 4*Beq;
Kvco = K/Kpd;
% ---------------------

LoopFilter = tf([tau2,1],[tau1,1], SamplesPerSecond);
VCOFilter = tf([Kvco],[1,0], SamplesPerSecond);

q_ki = zeros(NumberOfSymbols,1);
q_kq = zeros(NumberOfSymbols,1);

EsitmatedPhase = 0;
for i = 1:NumberOfSymbols
    % i at each step
    k = i-1;
    lower_bound = k*SamplesPerSymbol+1;
    upper_bound = (k+1)*SamplesPerSymbol;
    TimeInterval = t(1,lower_bound : upper_bound);
    
    [s1a,s2a] = PhaseDetector(s_r, EsitmatedPhase, omega_c, SamplesPerSecond, TimeInterval);
    
    s1b = s1a .* q_kq_prev;
    s2b = s2a .* q_ki_prev;

    % Phase error
    epsilon = s2b - s1b;

    LoopFilterOutput = lsim(LoopFilter, epsilon);
    NextEstimatedPhases = lsim(VCOFilter, LoopFilterOutput);

    % Demodulate the current symbol
    [qk_i, qk_q] = matched_demodulate( s_r , A_0, A_r, omega_c, const_phase_offset, NextEstimatedPhases, Ts, ...
                                       SamplesPerSecond);
    [ak_i, ak_q] = MLLDecision([qk_i, qk_q]', Symbols);
    
    %DecodedBits = decoder(,Symbols,SymbolBitVector);


end

end

