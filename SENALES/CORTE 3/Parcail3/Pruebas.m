%TEST INICIAL
clc
clear all

%Llamamos el archivo
file = "3notas.mp3";

%Sacamos la informacion de la cancion 
info = audioinfo(file)
[nota,Fs] = audioread(file);

nota = (nota(:,1)+nota(:,2))/2;

%Conseguir tiempo de la cancion
N = length(nota)
t = 0:seconds(1/Fs):seconds(info.Duration);
t = t(1:end-1);

%TRANSFORMADA
noteFFT = fft(nota) 
f = (0:N-1)*(Fs/2*N); 
noteFFT = noteFFT(1:floor(N/2));
f = f(1:floor(N/2))/10;
magnitude = abs(noteFFT); 

%%
%TEST INICIAL
clc;
clear all;
close all;

%Llamamos el archivo
file = "Parcail3/3notas_noise.mp3";

%Sacamos la informacion de la cancion 
info = audioinfo(file)
[nota,Fs] = audioread(file);

%nota = (nota(:,1)+nota(:,2))/2;
N = size(nota,1); % Determine total number of samples in audio file
df = Fs / N;
w = (-(N/2):(N/2)-1)*df;
y = fft(nota(:,1), N) / N; % For normalizing, but not needed for our analysis
y2 = fftshift(y);


n = 7;
beginFreq = 1000 / (Fs/2);
endFreq = 12000 / (Fs/2);
[b,a] = butter(n, [beginFreq, endFreq], 'bandpass');
noteOut = filter(b, a, nota);

subplot(3,1,1);
plot(1:N, nota);

subplot(3,1,2);
plot(w,abs(y2));

subplot(3,1,3);
plot(1:N, noteOut);

%%

alpha = 7;
cut = (523.25 - 0.000001)*10;
cutoff = (523.25 + 0.000001)*10;

nexttile
y = bandpass(nota,[500 501],Fs);

%TRANSFORMADA FILRTO
yFFT = fft(y) 
f = (0:N-1)*(Fs/N); 
yFFT = yFFT(1:floor(N/2));
f = f(1:floor(N/2))/10;
magnitudeY = abs(yFFT); 


nexttile
plot(t,nota)
xlabel('Time')
ylabel('Note')

nexttile
plot(0:0.5/length(magnitude):0.5-0.5/length(magnitude),magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Audio');

nexttile
plot(t,y)
xlabel('Time')
ylabel('Note')

nexttile
plot(f,magnitudeY)
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Audio FILTERED');


%% 
pitch(notaNoM,Fs)
notaFFT = fft(nota);
[maxValue, indexOfMaxValue] = max(notaFFT)


L = numel(t);                                   % Length Of Time & Signal Vectors
NotaFFT = fftshift(fft(nota)/L);
notaFFT = fft(nota);
[maxValue, indexOfMaxValue] = max(notaFFT)
stem(notaFFT)


