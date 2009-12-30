function [ s_di, s_dq , t] = transmit_filter( ak, bk, Ts, SamplesPerSecond)
%transmit_filter

NumberOfSymbols = length(ak);
dt = 1 / SamplesPerSecond;

% Pre-allocate arrays for better performance
% using round to prevent matlab warning about non-integers (although the
% result must be an integer)
s_di = zeros(1,round(Ts*NumberOfSymbols*SamplesPerSecond));
s_dq = zeros(1,round(Ts*NumberOfSymbols*SamplesPerSecond));
t = 0:dt:Ts*NumberOfSymbols; % calculate time sample values
t(:,end) = []; % Erase last time value (belongs to next symbol)

% Modulation
for i = 1:NumberOfSymbols
    a = ak(i);
    b = bk(i);
    
    % The bounds are really always integers, but floating point errors 
    % cause matlab to print warnings
    lower_bound = floor((i-1)*Ts*SamplesPerSecond+1);
    upper_bound = floor(i*Ts*SamplesPerSecond);
    
    TimeInterval = t(1,lower_bound:upper_bound);
    s_di_interval = a*ones(length(TimeInterval),1);
    s_dq_interval = b*ones(length(TimeInterval),1);
    % Write the samples into big array
    s_di(1,lower_bound:upper_bound) = s_di_interval;
    s_dq(1,lower_bound:upper_bound) = s_dq_interval;
end

end
