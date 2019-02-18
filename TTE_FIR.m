function [C,c,c_lp,Rl1,Rd1,L_tx,C_tx] = TTE_FIR(mi,sigma,r,Fs,fc,channelRessonance)

% The goal is to design a FIR filter that accurately modelates the TTE
% channel. We will use the fdesign.arbmagnphase function.
%
% n = filter order
% f = normalized frequency vector (0 to 1). Normalized by Fs/2?
% h = complex frequency response values

% design methods:'freqsamp','equiripple','firls'.

% -- Ressonant circuit parameters
Rl1 = 0;
Rd1 = 0;
L_tx = 0;
C_tx = 0;

N0 = 16e3;                                                                  % number of samples (from 0 to Fs)
N0_2 = N0/2;                                                                % number of samples (from 0 to Fs/2);

f = (Fs/(N0))*(0:1:(N0_2-1));                                               % frequency vector

if strcmp(channelRessonance,'classic')
    % --- Generating the Transfer Function
    delta = sqrt(2)./sqrt(2*pi*f*mi*sigma);
    T = r./delta;
    Fr = ((T.^2)./(pi*sigma*r^5)).*exp(-T).*sqrt(1 + 2*T + 2*T.^2).*exp(1i*(atan(T./(T+1)) - T - pi/2));
elseif strcmp(channelRessonance,'Ressonant')    
    [Fr,Rl1,Rd1,L_tx,C_tx] = VL_com_ressonancia( mi,sigma,r,Fs);
elseif strcmp(channelRessonance,'nonRessonant')
    [Fr,Rl1,Rd1,L_tx,C_tx] = VL_sem_ressonancia( mi,sigma,r,Fs );
end

Fr = Fr/max(abs(Fr));                                                      % Normalization

% OBS: aumentar a ordem do filtro de acordo com Fs (i.e. quanto maior Fs maior tem de ser a ordem do filtro N)!
N = 8*10;                                                                   % Filter order (max 8*60 min 8*10, always 8*n, n integer) 

d = fdesign.arbmagnphase('n,f,h',N,f/max(f),Fr); 

C = design(d,'freqsamp');

% --- Low pass equivalent
c = C.numerator;
t = (1/Fs)*(0:1:length(c) - 1);
c_lp = hilbert(c).*exp(-1i*2*pi*fc*t);

end