classdef Nota < handle
   properties 
      fileName
      signal
      duration
      sampleFrec
      frecFund
      octave
      note
   end
     
   methods 
       function object = setNote(object, file) 
            %SET FILENAME
            object.fileName = file;

            %SET SIGNAL AND SAMPLE FREC
            [tempSignal,tempSf] = audioread(file);

            if (tempSf == 48000)
                object.sampleFrec = tempSf;
                object.signal = tempSignal;
            else 
                fn = strrep(file,'.mp3','.wav');
                audiowrite(fn,tempSignal,48000);
                [object.signal,object.sampleFrec] = audioread(fn);
                delete(file);
                file = fn;
            end 
            
            %SET DURATION AND NUMBER OF CHANNELS
            info = audioinfo(file);
            object.duration = info.Duration;

            %CHECK IF MONO OR STEREO
            n = info.NumChannels;

            if(n > 1)
                sum = 0;
                for i = 1 : 1 : n
                        sum = sum + object.signal(:,i);
                end
                object.signal = sum/n;
            end 
       end
   end
end
