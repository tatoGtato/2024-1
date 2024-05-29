clc
clear all
close all

%%
lowCutoff = 300;  
highCutoff = 3400;
fs = 48000;

% Cargar el archivo de audio

%file = "baseDeDatosReChidori/6Voz_henrry.wav";
%[signal, Fs] = audioread(file);

recorder = audiorecorder(fs, 24, 1);
recordblocking(recorder, 3);
signal = getaudiodata(recorder, 'double');

signal = (signal ./ max(abs(signal), [], 1));
%%


[b, a] = butter(5, [lowCutoff, highCutoff] / (fs / 2), 'bandpass');
filteredSignal = filtfilt(b, a, signal);
sound(filteredSignal, fs);

%%
audiowrite("filti.wav",filteredSignal,48000);
%%
segmentlen = 100;
noverlap = 90;
NFFT = 128;

[S, F, T, P] = spectrogram(filteredSignal,segmentlen,noverlap,NFFT,fs,'yaxis');

% subplot(2,1,1)
% t = (0:length(filteredSignal)-1)/fs; 
% plot(t, signal);
% 
% subplot(2,1,2)
% t = (0:length(filteredSignal)-1)/fs; 
% plot(t, filteredSignal);

subplot(3,1,2)
spectrogram(filteredSignal,segmentlen,noverlap,NFFT,fs,'yaxis');

% subplot(4,1,3)
% spectrogram(filteredSignal,1024,512,1024,fs,'yaxis')


%%
P_dB = 10 * log10(P);

% Define a power threshold to identify high power regions (yellow regions)
powerThreshold = -50; % Adjust this value based on visual inspection

% Create a binary mask for regions above the power threshold
mask = P_dB > powerThreshold;

reconstructedSignal = zeros(size(filteredSignal));

for k = 1:length(T)
    % Apply the mask to the spectrogram column
    maskedSpectrum = S(:, k) .* mask(:, k);
    
    % Inverse FFT to convert back to time domain
    segment = ifft(maskedSpectrum, NFFT);
    
    % Determine the start and end indices for overlap-add
    startIdx = (k-1) * (segmentlen - noverlap) + 1;
    endIdx = startIdx + segmentlen - 1;
    
    % Ensure indices are within bounds
    if endIdx <= length(reconstructedSignal)
        % Overlap-add the segment into the reconstructed signal
        reconstructedSignal(startIdx:endIdx) = reconstructedSignal(startIdx:endIdx) + real(segment(1:segmentlen));
    end
end

% Normalize the reconstructed signal
reconstructedSignal = reconstructedSignal / max(abs(reconstructedSignal));

% Save or play the reconstructed signal
sound(reconstructedSignal, fs);

% Plot the original and reconstructed signals
t = (0:length(signal)-1)/fs;
figure;
subplot(2,1,1);
plot(t, signal);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, reconstructedSignal);
title('Reconstructed Signal');
xlabel('Time (s)');
ylabel('Amplitude');

audiowrite("filtiTHRESH.wav",reconstructedSignal,48000);

%%

startId = 1;
endId = 1;

sect = ones(size(reconstructedSignal));


for i = 1: length(reconstructedSignal)
    if(reconstructedSignal(i) ~= 0)
        startId = i;
        break
    end
end 
    
for i = startId: length(reconstructedSignal)
    if(reconstructedSignal(i:i+1000) == 0)
            reconstructedSignal(i:i+1000)
            endId = i;
            break
    end
end 

t(startId)
t(endId)
croppedT = t(startId : endId);
croppedSig = reconstructedSignal(startId : endId);

subplot(2,1,1);
hold on
plot(t, reconstructedSignal);
area(t(t> t(startId) &t< t(endId) ),sect(t>t(startId)&t<t(endId)),"FaceColor","r","FaceAlpha",0.5)
area(t(t> t(startId) &t< t(endId) ),-sect(t>t(startId)&t<t(endId)),"FaceColor","r","FaceAlpha",0.5)
hold off

subplot(2,1,2);
plot(croppedT, croppedSig)
%%
t1 = t(startId)
t2 = t(endId)

dt = 1/fs;
I0 = round(t1/dt);
Iend = round(t2/dt);
x = reconstructedSignal(I0:Iend);

t = 0:dt:length(x)*dt-dt;

plot(t,x)

c = cceps(x);

[C,I] = max(c)


plot(t,c)
xlabel('ms')

hold on
plot(t(I),c(I),'o')
hold off
