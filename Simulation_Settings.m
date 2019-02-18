% -- Simulation Settings

% ----------------------Control Parameters---------------------------------
channelCoding = 0;                                                          % enables the use of Reed Solomon channel coding (0 or 1)
modulationScheme = 'PSK';                                                   % PSK,PAM
pulseShaping = 'RCOS';                                                      % RCOS,modDuobin
noiseType = 'AWGN';                                                         % AWGN,TTE,none
channelType = 'TTE';                                                        % ideal, TTE
channelRessonance = 'classic';                                              % classic, Ressonant, nonRessonant (obs: classic = nonRessonant)
equalization = 1;                                                           % enables MLSE equalization


if strcmp(modulationScheme,'PAM')
    baseBand = 1;
elseif strcmp(modulationScheme,'PSK')    
    baseBand = 0;
end
if strcmp(channelType,'ideal')
    equalization = 0;                                                       % If the channel is ideal, no need for equalization
end
% ----------------------Control Parameters---------------------------------

% ----------------Modulation Parameters and objects------------------------
modOrder = 2;
pamMod = comm.PAMModulator(modOrder,'BitInput',true);
pamDemod = comm.PAMDemodulator(modOrder,'BitOutput',true);
pskMod = comm.PSKModulator(modOrder, 'BitInput',true,'PhaseOffset',0);
pskDemod = comm.PSKDemodulator(modOrder, 'BitOutput',true,...
            'DecisionMethod','Hard decision','PhaseOffset',0);

% --- Choosing which mod/demod to use
if strcmp(modulationScheme,'PAM')

    hModulator = pamMod;
    hDemod = pamDemod;

elseif strcmp(modulationScheme,'PSK')
    hModulator = pskMod;
    hDemod = pskDemod;
end
% ----------------Modulation Parameters and objects------------------------

% -------------------Reed Solomon channel coding---------------------------
m_rs = 8;                                                                   % Defines the Galois Field in which the code operates
t_rs = 102;                                                                 % Number of symbol errors the code can correct
N_rs = 2^m_rs - 1;                                                          % Length of the output (codeword)block
K_rs = N_rs - 2*t_rs;                                                       % Length of the input (message)block
R_rs = K_rs/N_rs;                                                           % Rate of the error-correcting code

if channelCoding
    enc = comm.RSEncoder(N_rs,K_rs,'BitInput',1);                           % Encoder object
    dec = comm.RSDecoder(N_rs,K_rs,'BitInput',1);                           % Decoder object
else
    R_rs = 1;
end
% -------------------Reed Solomon channel coding---------------------------

% -------------------------General parameters------------------------------
Numloops = 1e0;                                                             % Number of simulation loops
NumDataBits = 1e6;                                                          % Number of Data Bits to be simulated
numBlocks = ceil(NumDataBits/(m_rs*K_rs*log2(modOrder)));                   % Number of simulated code blocks
numBits = numBlocks*m_rs*K_rs;                                              % Number of simulated message bits.
numCodedBits = numBlocks*m_rs*N_rs;                                         % Number of coded bits
fb = 75000;                                                                 % Information bit rate
Rs = fb/(R_rs*log2(modOrder));                                              % Symbol rate (in bauds)
EbN0 = [0 1 2 3 4 5 6 7 8 9 10];                                            % Eb/N0 vector
hError = comm.ErrorRate;                                                    % Error object

% -------------------------General parameters------------------------------

% -----------Pulse Shaping filters and Sampling parameters-----------------
Nsym = 24;                                                                  % Duration of the filter impulse response in number of symbols
sampsPerSym = 8;                                                            % Oversampling factor
beta = 0.5;                                                                 % RCOS Roll - off.
B_base = Rs*(1+beta)/2;                                                     % baseband bandwidth
B_mod = Rs*(1+beta);                                                        % passband bandwidth
if baseBand
    BW = B_base;                                                            % Signal Bandwidth
else
    BW = B_mod;
end

if Rs < 3750
    Fs = 3750*sampsPerSym;                                                  % Por que fazer isso mesmo?
else
    Fs = Rs * sampsPerSym;                                                  % Sampling frequency (samples/s).
end

[filterTx,filterRx,fltDelay] = RCOS_pulse_setup(Nsym,sampsPerSym,Fs,beta ); % Filter RCOS

% --- Duobin filters
F_duobin = (1/sqrt(2))*[1; 0; -1];
%F_duobin = (1/sqrt(2))*[1; 1];
C_duobin = dfilt.dffir(upsample(F_duobin,sampsPerSym));
c_mod_duobin = upsample(F_duobin,sampsPerSym);

% -----------Pulse Shaping filters and Sampling parameters-----------------

% --------------------Channel/Noise parameters-----------------------------
mi = 4*pi*1e-7;                                                             % magnetic permeability (in H/m)
sigma = 0.01;                                                               % conductivity (in S/m)
r = 200;                                                                    % distance between transmitter and receiver (in m)
f_opt = 1*(16/(2*pi*sigma*mi*r^2));

if B_mod < 2*f_opt
    fc = f_opt;
else
    fc = 1*(1+beta)*Rs/2;
end
Rs_max = 2*fc/(1 + beta);                                                    % Maximum symbol rate (in bauds).

% --- TTE Channel (FIR system)   
[C,c,c_lp,Rl1,Rd1,L_tx,C_tx] = TTE_FIR(mi,sigma,r,Fs,fc,channelRessonance);
if baseBand
    C_TTE = c;                                                              % normal channel
else
    C_TTE = c_lp;                                                           % low pass equivalent
end

% --- VLF Noise model(Field-Lewinstein)
imp_2 = 7;                                                                  % square of the impulsivity
a = 0.7;                                                                    % Weibull shape parameter 
% --------------------Channel/Noise parameters-----------------------------

%--------------------------Equalization------------------------------------

% White Matched Filter for TTE channel
% (use when baseband = 1 && RCOS)
B = conv(filterTx.numerator,c);
B_inv = conv(fliplr(c),filterRx.numerator);
D = conv(B,B_inv);
E = downsample(D,sampsPerSym);
F = firminphase(E);

% White Matched Filter for low pass equivalent TTE channel 
%(use when baseband = 0 && RCOS)
B_lp = conv(filterTx.numerator,c_lp);
B_inv_lp = conv(fliplr(conj(c_lp)),filterRx.numerator);
D_lp = conv(B_lp,B_inv_lp);
E_lp = downsample(D_lp,sampsPerSym);
F_lp = complexfirminphase(E_lp);
 
% White Matched Filter for modDuobin pulse + TTE channel    
% (usar quando baseband = 1 && modDuobin)

%OBS: Bugzinho com F_mod_duobin. Para algumas taxas, WMF não é calculado
%corretamente, por causa de zeros sobre o círculo unitário. Verificar qual 
% é o filtro correto: F_mod_duobin ou F2_mod_duobin

A_mod_duobin = conv(filterTx.numerator,C_duobin.numerator);
B_mod_duobin = conv(A_mod_duobin,c);
B_inv_mod_duobin = conv(fliplr(c),conv(filterRx.numerator,fliplr(C_duobin.numerator)));
D_mod_duobin = conv(B_mod_duobin,B_inv_mod_duobin);
E_mod_duobin = downsample([0 D_mod_duobin],sampsPerSym);
E_mod_duobin = E_mod_duobin(2:end);
F_mod_duobin = firminphase(E_mod_duobin);
F2_mod_duobin = firminphase(E_mod_duobin,2);

if Rs == 20000
    F_mod_duobin = F2_mod_duobin;
end
    
%------------------------Equalization--------------------------------------

% ---------------------Theoretical performances----------------------------
EbN0_abs = 10.^(EbN0./10);                                                  
bpsk_teorico = qfunc(sqrt(2.*EbN0_abs));
duobinario = 2*(3/4)*qfunc(sqrt(((pi/4)^2)*2*EbN0_abs));
% ---------------------Theoretical performances----------------------------
