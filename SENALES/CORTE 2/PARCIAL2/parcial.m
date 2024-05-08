%%
%TEST INICIAL
clc
clear all

%Llamamos el archivo
file = "LalaTrim45.mp3";

%Sacamos la informacion de la cancion 
info = audioinfo(file)
[cancion,Fs] = audioread(file);

%Sacamos frecuencia de muestreo de la cancion
frecMuestreo = info.SampleRate
muestrasTotales = info.TotalSamples
duracion = info.Duration

%Como la cancion es estereo pasarla a mono
sum = 0;

for i = 1 : 1 : info.NumChannels
    sum = sum + cancion(:,i);
end

cancionforMono = sum/info.NumChannels;
cancionMono = (cancion(:,1)+cancion(:,2))/2;

if(cancionMono == cancionforMono)
    "true"
end

%Conseguir tiempo de la cancion
t = 0:seconds(1/Fs):seconds(info.Duration);
t = t(1:end-1);

%Para plotear 
tiledlayout(2,1)
% Top plot
nexttile
plot(t,cancion)
xlabel('Time')
ylabel('Multi channel Audio Signal')

% Bottom plot
nexttile
plot(t,cancionMono)
xlabel('Time')
ylabel('Single channel Audio Signal')

%%
%Test con dos canciones 

clc
clear all

%Importar y tomar informacion de las canciones y el imput
song1 = "superBaseDeDatosSuperSeguraConSuperDatos/Lala.mp3";
info = audioinfo(song1)
[cancion1,Fs1] = audioread(song1);

song2 = "superBaseDeDatosSuperSeguraConSuperDatos/takeTrain.mp3";
info2 = audioinfo(song2)
[cancion2,Fs2] = audioread(song2);

input = "LalaTrim45.mp3";
infoInp = audioinfo(input)
[cancionI,FsI] = audioread(input);

%%
a.a = cancion1

a.b = cancion2



%%

%Ponerlas mono
cancion1Mono = (cancion1(:,1)+cancion1(:,2))/2;
cancion2Mono = (cancion2(:,1)+cancion2(:,2))/2;
inputMono = (cancionI(:,1)+cancionI(:,2))/2;

%Graficar las tres pistas de sonido

t1 = 0:seconds(1/Fs1):seconds(info.Duration);
t1 = t1(1:end-1);

t2 = 0:seconds(1/Fs2):seconds(info2.Duration);
t2 = t2(1:end-1);

t3 = 0:seconds(1/FsI):seconds(infoInp.Duration);
t3 = t3(1:end-1);


tiledlayout(3,1)
% Song1 plot
nexttile
plot(t1,cancion1Mono)
xlabel('Time')
ylabel('Melody 3')
grid on
xlim(seconds([0 180]))
ylim([-2 2])

% Song2 plot
nexttile
plot(t2,cancion2Mono)
xlabel('Time')
ylabel('Take the A train')
grid on
xlim(seconds([0 180]))
ylim([-2 2])

% Song3 plot
nexttile
plot(t3,inputMono)
xlabel('Time')
ylabel('Melody 3 fragment')
grid on
xlim(seconds([0 180]))
ylim([-2 2])


%%
clear all
clc

[raton,FsR] = audioread("aventurasRaton.mp3");
[cancion2,Fs2] = audioread("Melody3.mp3");

FsR
Fs2

%audiowrite('aventurasRaton48.wav',raton,Fs2);

[cancion2,F] = audioread("aventurasRaton48.wav");
F
%%
%Sacar la correlacion 

[C1,lag1] = xcorr(cancion1Mono,inputMono);        
[C2,lag2] = xcorr(cancion2Mono,inputMono);        

figure

ax(1) = subplot(2,1,1); 
plot(lag1/Fs1,C1,"k")
ylabel("Amplitude")
grid on
title("Cross-Correlation entre Lala (Cancion coincidencia) y el input")

ax(2) = subplot(2,1,2); 
plot(lag2/Fs2,C2,"r")
ylabel("Amplitude") 
grid on
title("Cross-Correlation entre Take the A train y el input")
xlabel("Time(s)") 
axis(ax(1:2),[0 200 -700 700])

max(abs(C1))
max(abs(C2))

%%
%Encontrar el delay
%En este punto ya sabemos que la cancion1 es la misma cancion del input

[~,I] = max(abs(C1))
SampleDiff = lag1(I);
timeDiff = SampleDiff/Fs1

%%
%Matchear las funciones
%Ya sabemos el delay solo falta matchearlas

tiledlayout(2,1)
% Song1 plot
nexttile
plot(t1,cancion1Mono)
xlabel('Time')
ylabel('Melody 3')
grid on
xlim(seconds([0 160]))
ylim([-2 2])

tf = 0+seconds(timeDiff):seconds(1/FsI):seconds(infoInp.Duration)+seconds(timeDiff);
tf= tf(1:end-1);

% Song3 plot
nexttile
plot(tf,inputMono)
xlabel('Time')
ylabel('Melody 3 fragment')
grid on
xlim(seconds([0 160]))
ylim([-2 2])

