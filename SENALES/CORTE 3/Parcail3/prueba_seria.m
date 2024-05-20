%%PRIMERO QUITAR EL RUIDO DE FONDO

clc
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

N = size(noteOut,1); % Determine total number of samples in audio file
df = Fs / N;
w2 = (-(N/2):(N/2)-1)*df;
y = fft(noteOut(:,1), N) / N; % For normalizing, but not needed for our analysis
y3 = fftshift(y);


subplot(4,1,1);
plot(1:N, nota);

subplot(4,1,2);
plot(w,abs(y2));

subplot(4,1,3);
plot(1:N, noteOut);

subplot(4,1,4);
plot(w2,abs(y3));

%% DEJAR SOLO LAS NOTAS usando el pitch 
clc
close all;

filteredPitch = pitch(noteOut,Fs);

length(filteredPitch)
length(noteOut)

constant_indices = [];
tolerance = 1e-6; 

for k = 1:length(filteredPitch) - 6
    if all(abs(diff(filteredPitch(k:k+6))) < tolerance)
        constant_indices = [constant_indices, k];
    end
end

onlyNotesPitch = filteredPitch(constant_indices(1):constant_indices(end));
d = constant_indices(1):constant_indices(end);


subplot(2,1,1);
pitch(noteOut,Fs)

subplot(2,1,2)
plot(d, onlyNotesPitch);


%% SEPARAR NOTA POR NOTA
clc 

uA = unique(onlyNotesPitch)


