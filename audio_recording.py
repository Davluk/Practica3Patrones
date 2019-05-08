import sounddevice as sd
import numpy as np
import scipy.io.wavfile as wav
import csv
import os

file = 'diez'
fs=16000
duration = 2  # seconds

if(not os.path.isdir('recs/'+file)):
        os.mkdir('recs/'+file)
for index in range(15):
    myrecording = sd.rec(duration * fs, samplerate=fs, channels=2)
    print (index," Recording Audio")
    sd.wait()
    print ("Audio recording complete , Play Audio")
    sd.play(myrecording, fs)
    sd.wait()
    print ("Play Audio Complete")
    with open('recs/'+file+'/'+file+str(index)+".csv", "w") as write_file:
        writer = csv.writer(write_file)
        writer.writerows(myrecording)


