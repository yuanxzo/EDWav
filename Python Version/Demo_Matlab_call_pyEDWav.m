% This is a simple demonstration program that calls EDWav.py on the Matlab platform
% Require your MATLAB to load Python modules correctly. Use >> pyenv for testing.

clear

% load EDWav.py library 
EDWav = py.importlib.import_module('EDWav');  

% load a Matlab formatted test seismic waveform data
load('Demo_data/Demo_data.mat');

% evaluate this waveform
obj=EDWav.EDWav.evaluate(wvfm,FS=Fs);

% Convert numpy.ndarray data format to double format
obj1.A = double(obj.A);
obj1.B = double(obj.B);
obj1.C = double(obj.C);
obj1.freqs = double(obj.freqs);
obj1.proxy = double(obj.proxy);


figure;
imagesc(obj1.freqs,obj1.freqs,obj1.B);
