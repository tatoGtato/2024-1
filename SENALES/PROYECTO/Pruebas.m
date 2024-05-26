clc
clear all
close all

%% Se crea un objeto recorder, para grabar (Pero por ahora solo vamos a )

% fs = 48000;
% recorder = audiorecorder(fs, 24, 1);
% 
% disp("Porfavor graba tu voz: ");
% drawnow();
% recordblocking(recorder, 5);
% data = getaudiodata(recorder, 'double');

file = "voiceTest.mp3";
info = audioinfo(file)
[data,fs] = audioread(file);

t = (0:length(data)-1)/fs; 
plot(t, data);
xlabel('Time (s)');
ylabel('Amplitude');
title('Time-domain Signal');

%% Con lo encontrado en la FTT vamos a extaer las features de la voz y ponerlas en una base de datos
f = getFeaturesVoz(data);

persona = input("Pon el nombre de la persona: ");
try 
    load baseDeDatosReChidoriDB
    if exist('frec', 'var') && exist('name', 'var')
        frec = [frec; f];
        name = [name; persona];
    else
        frec = f;
        name = persona;
    end
    save baseDeDatosReChidoriDB frec name
catch 
    frec = f;
    name = persona;
    save baseDeDatosReChidoriDB frec name
end 
msgbox("La nota de voz se ha registrado")

%% TEST DEL MODELO (Clasificacion) encontrar las features que ams separezcan a las base de datos

f = getFeaturesVoz(data)

load baseDeDatosReChidoriDB
D = [];
for i = 1: size(frec,1)
    d = sum(abs(frec(i) - f));
    D = [D d];
end 

sm = inf;
ind = -1;
for i = 1 : length(D)
    if (D(i) < sm)
        sm = D(i);
        ind = i;
    end 
end 
personaDetectada = name(ind);
disp("La persona detectada es la: ");
personaDetectada


%% FUnciones

function xPitch = getFeaturesVoz(data)
    f = real(fft(data));
    plot(abs(f));
    maxFrec = max(f);
    xPitch = find(real(f) == maxFrec,1);
end 




