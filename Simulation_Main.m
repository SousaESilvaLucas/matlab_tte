%% Through - the - earth link simulation
clear
clc
tic
% --- Run Simulation_Settings.m
Simulation_Settings 

ber = zeros(length(EbN0),1);                                                % Initialization of BER vector
for indice = 1:length(EbN0)
     hError.reset                                                           % Initialization of error counter object
     for counter = 1:Numloops
        fprintf('Eb/N0 = %f \nLoop number %i\n',EbN0(indice),counter)
        toc
        % --- Randomly generating bits
        messageBits = randi([0 1],numBits,1); 

        if channelCoding
        % --- Encoding the generated bits
            codedBits = step(enc,messageBits);
        else
            codedBits = messageBits;
        end

        % --- Modulate the data
            modData = step(hModulator, codedBits);

        % --- Pulse shaping + Transmission through channel

        % --- Oversampling and Tx RCOS filtering
            yo = filter(filterTx, upsample([modData; zeros(Nsym/2,1)],...
                                                        sampsPerSym));
            yo = yo(fltDelay+1:end);                                        % Removing filter delay

        % --- Duobin Correlation filter
        if strcmp(pulseShaping,'modDuobin')     
           yt = filter(C_duobin,yo);
        elseif strcmp(pulseShaping,'RCOS')
            yt = yo;
        end

        % --- Signal passes through channel
        if strcmp(channelType,'ideal')
            y_channel = yt;
        elseif strcmp(channelType,'TTE')
            y_channel = filter(C_TTE,1,yt);
        end

        P_signal = 10*log10(mean(abs(y_channel).^2));                       % Signal energy
        % --- Calculates SNR
        if channelCoding
            SNR = EbN0(indice) + 10*log10(R_rs*log2(modOrder)) - 10*log10(sampsPerSym);           
        else
            SNR = EbN0(indice) + 10*log10(log2(modOrder)) - 10*log10(sampsPerSym);           
        end

        P_noise = P_signal - SNR;                                          % Noise power

        % --- Corruption with noise
        if strcmp(noiseType,'AWGN')
            if isreal(y_channel)
                y_noise = awgn(y_channel,SNR + 3,'measured');
            else
                y_noise = awgn(y_channel,SNR,'measured');
            end
        elseif strcmp(noiseType,'TTE')
            if isreal(y_channel)
                noise = real(noiseTTE(a,imp_2,length(y_channel),P_noise));
                y_noise = y_channel + noise;
            else
                noise = noiseTTE(a,imp_2,length(y_channel),P_noise);
                y_noise = y_channel + noise;
            end
        elseif strcmp(noiseType,'none')
            y_noise = y_channel;
        end

        % --- Matched Filtering (Channel)
        if strcmp(channelType,'ideal')
            y_receiver = y_noise;
        elseif strcmp(channelType,'TTE')    
            y_receiver = filter(conj(fliplr(C_TTE)),1,[y_noise; zeros(length(C_TTE)-1,1)]);
            y_receiver = y_receiver(length(C_TTE):end);                         % Removing filter delay
        end
        % --- Matched Filtering (RCOS)
        yr = filter(filterRx,[y_receiver; zeros(Nsym*sampsPerSym/2, 1)]);
        yr = yr(fltDelay+1:end);                                                % Removing filter delay

        % --- Matched Filtering (Duobin)
        if strcmp(pulseShaping,'modDuobin')
            yr = filter(fliplr(C_duobin.numerator),1,[yr; zeros(length(C_duobin.numerator)-1,1)]);
            yr = yr(length(C_duobin.numerator):end);                            % Removing filter delay
        end

        % --- Signal Detection
        if strcmp(modulationScheme,'PAM')
            if strcmp(pulseShaping,'modDuobin')
                if equalization
        % --- Noise Whitening Filter
                    receivedSamples = filter(1,F_mod_duobin,flipud(downsample(yr,sampsPerSym)));
                    receivedSamples = flipud(receivedSamples);
        % --- MLSE equalization
                    eqSamples = mlseeq(receivedSamples,F_mod_duobin(1:7),[1 -1],30,'rst');
                else
        % --- Noise Whitening Filter
                    receivedSamples = filter(1,F_duobin,flipud(downsample(yr,sampsPerSym)));
                    receivedSamples = flipud(receivedSamples);
        % --- MLS Detection (no channel equalization)
                    eqSamples = mlseeq(receivedSamples,F_duobin,[1 -1],20,'rst');
                end

            elseif strcmp(pulseShaping,'RCOS')
                if equalization
        % --- Noise Whitening Filter
                    receivedSamples = filter(1,F,conj(flipud(downsample(yr,sampsPerSym))));
                    receivedSamples = conj(flipud(receivedSamples));
        % --- MLSE equalization
                    eqSamples = mlseeq(receivedSamples,F(1:7),[1 -1],30,'rst');
                else
                    eqSamples = downsample(yr,sampsPerSym);
                end
            end
        elseif strcmp(modulationScheme,'PSK')
            if equalization
        % --- Noise Whitening Filter
                receivedSamples = filter(1,F_lp,conj(flipud(downsample(yr,sampsPerSym))));
                receivedSamples = conj(flipud(receivedSamples));   
        % --- MLSE equalization
                eqSamples = mlseeq(receivedSamples,F_lp(1:7),[1 -1],30,'rst');
            else
                eqSamples = downsample(yr,sampsPerSym);
            end
        end

        % --- Demodulation
        if isreal(eqSamples)
            receivedBits = step(hDemod,complex(eqSamples));
        else
           receivedBits = step(hDemod,eqSamples);
        end


        if channelCoding
        % --- Decoding
            decodedBits = step(dec,receivedBits);
        else
            decodedBits = receivedBits;
        end      

        % --- Bit Error Rate
        errorStats  = step(hError, messageBits, decodedBits);
        ber(indice) = errorStats(1);
     end
end
   
      
%% Gráficos

% --- BER plot
figure(2)
semilogy(EbN0,ber,'o-')
xlabel('Eb/N0 (dB)')
ylabel('BER')
title('BER em função de Eb/N0')
grid on
hold on
% plot(EbN0,duobinario,'-r')
plot(EbN0,bpsk_teorico,'-r')

% visualisação da banda do sinal
%     figure(7)
%     f = (0:1:length(y_bandpass)-1)*(Fs/length(y_bandpass));
%     plot(f,abs(Y_BANDPASS)/max(abs(Y_BANDPASS)))
%     hold on
%     plot(f,abs(H_channel),'-r')
%     xlabel('frequência (hz)')

toc