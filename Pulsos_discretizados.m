clc
clear
numBits = 25;
pulseShaping = 'PRS';                                                      % RCOS ou PRS (partial response signal)
Nsym = 24;                                                                  % Duration of the filter impulse response in number of symbols
sampsPerSym = 8;                                                            % Oversampling factor
beta = 0.5;                                                                 % RCOS Roll - off.
Rs = 2000;                                                                  % Symbol Rate
Fs = Rs * sampsPerSym;                                                      % Sampling Frequency
[filterTx,filterRx,fltDelay] = RCOS_pulse_setup(Nsym,sampsPerSym,Fs,beta ); % RCOS Filter
% -- Partial response filters 
F_prs = [1; 0; -1];                                                         % [1;0;-1] = modified duobinary; [1;1] = duobinary.                                              
C_prs = dfilt.dffir(upsample(F_prs,sampsPerSym));

% -- Symbols to be transmitted
messageBits = randi([0 1],numBits,1);                                       % random bits
messageBits = [1 0 0 1 0 1 1 0]';
modData = messageBits*2 - 1;


% --- Oversampling and Tx RCOS filtering
    yo = filter(filterTx, upsample([modData; zeros(Nsym/2,1)],...
                                                sampsPerSym));
    yo = yo(fltDelay+1:end);                                        % Removing filter delay

% --- Duobin Correlation filter
if strcmp(pulseShaping,'PRS')     
   yt = filter(C_prs,yo);
elseif strcmp(pulseShaping,'RCOS')
    yt = yo;
end

t = Fs*(0:1:length(yt)-1);
figure(99)
stem(t,yt)
