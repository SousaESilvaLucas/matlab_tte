function [ sqrtRcosFlt,sqrtRcosFltRcv,fltDelay] = RCOS_pulse_setup(Nsym,sampsPerSym,Fs,beta )

% 'Square Root Raised Cosine'filter design
shape = 'Square Root Raised Cosine';                                        % Pulse shape
%Fs = R * sampsPerSym;                                                       % Sampling frequency (samples/s).
fltDelay = (Nsym /2)*sampsPerSym;                                           % Filter delay
%BW = (1 + beta)*R;                                                           % Signal bandwidth

sqrtRcosSpec = fdesign.pulseshaping(sampsPerSym, shape, 'Nsym,beta',...
                                                            Nsym, beta,Fs);
sqrtRcosFlt = design(sqrtRcosSpec);                                         % Transmission filter.
sqrtRcosFltRcv = design(sqrtRcosSpec);                                      % Reception filter.

% Normalizando os filtros de transmissão e recepção.
normFact = max(sqrtRcosFlt.Numerator);
sqrtRcosFlt.Numerator = sqrtRcosFlt.Numerator/normFact;
sqrtRcosFltRcv.Numerator = sqrtRcosFltRcv.Numerator*(normFact*sampsPerSym);

% fvtool(sqrtRcosFlt)
% Filter order: Nsym*sampsPerSym + 1
% Filter delay: Nsym*sampsPerSym/2

end

