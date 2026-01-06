% This is a simple demonstration program that calls EDWav.py on the Matlab platform
% Require your MATLAB to load Python modules correctly. Use >> pyenv for testing.

clear

% load EDWav.py library 
EDWav = py.importlib.import_module('EDWav');  

% load a Matlab formatted test seismic waveform data
load('Demo_data/Demo_data.mat');

% If any code in EDWav.py has been modified, you can run this command to reload the module
% EDWav = py.importlib.reload(EDWav); 

% evaluate this waveform
obj=EDWav.EDWav.evaluate(wvfm,Fs);

% Convert numpy.ndarray data format to double format
obj1.A = double(obj.A);
obj1.B = double(obj.B);
obj1.C = double(obj.C);
obj1.freqs = double(obj.freqs);
obj1.proxy = double(obj.proxy);


figure;
imagesc(obj1.freqs,obj1.freqs,obj1.B);

% Here, we demonstrate the use of the 'EDWav.EDWav.evaluate' function on the MATLAB platform. 
% While the 'EDWav.EDWav.evaluate_sliding' function is not presented, it can be used in a similar manner.

