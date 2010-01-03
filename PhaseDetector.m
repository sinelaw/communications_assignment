function [ s1, s2 ] = PhaseDetector( s_r, est_phase, A_0, omega_c, constant_phase, t)
%PhaseDetector Detects the phase in the incoming signal s_r containing ONE symbol using the estimate phase est_phase and carrier
%frequency omega_c. A_0 is local oscillator amplitude. f_s is sampling frequency for simulation.

s1 =  A_0.* cos(omega_c*t + est_phase + constant_phase) .* s_r;
s2 = -A_0.* sin(omega_c*t + est_phase + constant_phase) .* s_r;


end

