function [ VL,Rl1,Rd1,L_tx,C_tx ] = VL_com_ressonancia( mu,sigma,r,Fs )
% r = distância entre as antenas
% sigma = condutividade do solo
% mi = permeabilidade magnética
% VL = tensão recuperada na antena receptora

% Considera-se que a corrente (I_Tx) que controlamos (sinal de entrada) = 1 A. Dessa forma VL = impedância de transferência do sistema
% Os parâmetros que modificarão o sistema com ressonância são os paramêtros
% de antena (Antenna parameters) e os valores dos resistores  de dumping e
% os capacitores do circuito (C_tx, C_rx, Rd1 e Rd2).


% -- Antenna parameters

raio_Tx=20;
N_Tx=1;
raio_Rx=1;
N_Rx=50;
raio_wire = 1e-3;
S_Tx=pi*raio_Tx^2;
S_Rx=pi*raio_Rx^2;

%------------------------------ LUCAS-------------------------------------%

% -- TTE Cannel construction
N0 = 16e3;                                                                  % number of samples (from 0 to Fs)
N0_2 = N0/2;                                                                % number of samples (from 0 to Fs/2);

f = (Fs/(N0))*(0:1:(N0_2-1));

dF = Fs/N0;
Fo = 0;
w=2*pi*f;

    % --- Generating the Transfer Function
    delta = sqrt(2)./sqrt(2*pi*f*mu*sigma);
    T = r./delta;
    Fmic = ((T.^2)./(pi*sigma*r^5)).*exp(-T).*sqrt(1 + 2*T + 2*T.^2).*exp(1i*(atan(T./(T+1)) - T - pi/2));
    %%% Canal  - medidas dos anos 70
    %Fmic = empirical_model_70(r,w,sigma,mu);  % Não é usado!

%------------------------------ LUCAS-------------------------------------%


%f_opt = 2.8312^2/(pi*mu*sig*r^2);
%f_opt = 1e3;
[~,s2] = max(Fmic); 
f_opt = f(s2);

%%%%%

% Circuit paramesters
Vs=1.363;%Vs=50; 
Rs=50; %Rs=1e-10; %50;
RL=1e6; 
sig_copp=5.85e7;
Rd1=0;Rd1=10; %[3 original]Tx damping load (amortecimento)
Rd2=0;Rd2=10; % [260 original]Rx damping load (amortecimento)

Rl1=2*N_Tx*raio_Tx/(sig_copp*raio_wire^2);   %Rl1=0.001;
Rl2=2*N_Rx*raio_Rx/(sig_copp*raio_wire^2);   %Rl1=0.001;

L_tx=N_Tx^2*mu*raio_Tx*log(8*raio_Tx/raio_wire -2) ; %L_tx=1.15e-6; 
L_rx=N_Rx^2*mu*raio_Rx*log(8*raio_Rx/raio_wire -2) ; %L_rx=1.15e-6; 

Cs = 1/(L_tx*(2*pi*f_opt)^2); Cs=5e+30;      %Series capacitor in Tx % For current source C-Series in Tx is USELESS 
C_tx = 1/(L_tx*(2*pi*f_opt)^2); %C_tx=1e-30;   %parallel capacitor in Tx
C_rx = 1/(L_rx*(2*pi*f_opt)^2); %C_rx=1e-30;   %parallel capacitor in Rx

%%% Circuito Tx - VOLTAGE SOURCE
K =1+1i*w*(Rl1+Rd1)*C_tx-w.^2*L_tx*C_tx; 
Z_1 = (Rl1+Rd1) + 1i*w*L_tx + 1./(1i*w*Cs); Z_2 = 1./(1i*w*C_tx);
Z_volt = Z_1.*Z_2./(Z_1+Z_2);
V=Vs*(1i*w*(Rl1+Rd1)*Cs-w.^2*L_tx*Cs)./(K+1i*w*Cs.*(K*Rs+(Rl1+Rd1))-w.^2*L_tx*Cs); 
%%%%%%%%%%%%%%%%%%%

%%% Circuito Tx - CURRENT SOURCE
Z_L = (Rl1+Rd1)+1i*w*L_tx; 
Z_C = 1./(1i*w*C_tx);
Z_cur2 = Z_L.*Z_C./(Z_L+Z_C); 
Z_cur= (Rl1 + Rd1 + 1i*w*L_tx)./(1i*w*Rl1*C_tx + 1i*w*Rd1*C_tx - C_tx*L_tx*w.^2 + 1);
%I_Tx = Vs/Rs;

I_Tx = 1;                       % LUCAS - Coloca-se I_Tx = 1 para encontrar a impedância de transferência
V=Z_cur*I_Tx; %coment this line if VOLTAGE SOURCE
%%%%%%%%%%%%%%%%%%%

%%% Tx Current
I_L=V./((Rl1+Rd1)+1i*w*L_tx); 
%%%%%%%%%%%%%%%%%%%

Zmic = N_Tx*N_Rx*S_Tx*S_Rx*Fmic;
Md = N_Tx*S_Tx*I_L;

%%% Circuito Rx
Vl = Zmic .* I_L; %Vl=1;
VL = Vl*RL./((RL+(Rl2+Rd2)) + 1i*w*(L_rx+C_rx*(Rl2+Rd2)*RL) - w.^2*L_rx*C_rx*RL); 

%Calculos de largura de banda BW e fator de qualidade Q
%Tx

%spec = abs(Md);
spec = abs(I_L);


[BW_ind, peak_ind, peak_val_Tx] = BandWidth(spec,1);
BW_Tx = dF*BW_ind; Fpeak_Tx = (peak_ind-1)*dF + Fo; Q_Tx = Fpeak_Tx/BW_Tx; 
%Canal
spec = abs(Fmic);
[BW_ind, peak_ind, peak_val_canal] = BandWidth(spec,1);
BW_canal = dF*BW_ind; Fpeak_canal = (peak_ind-1)*dF + Fo; Q_canal = Fpeak_canal/BW_canal; 
%Rx
VL_fonte_flat = RL./((RL+(Rl2+Rd2)) + 1i*w*(L_rx+C_rx*(Rl2+Rd2)*RL) - w.^2*L_rx*C_rx*RL) *N_Rx*S_Rx;
spec = abs(VL_fonte_flat);
[BW_ind, peak_ind, peak_val_Rx] = BandWidth(spec,1);
BW_Rx = dF*BW_ind; Fpeak_Rx = (peak_ind-1)*dF + Fo; Q_Rx = Fpeak_Rx/BW_Rx; 
%Composite system
spec = abs(VL);
[BW_ind, peak_ind, peak_val, left_ind, right_ind] = BandWidth(spec,1);
BW = dF*BW_ind; Fpeak = (peak_ind-1)*dF + Fo; Q = Fpeak/BW;

Pow_Dens = (Rl1+Rd1)*abs(I_L).^2; 
Pow = sum(Pow_Dens(left_ind:right_ind))*dF;

end

