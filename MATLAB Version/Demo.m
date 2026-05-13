% This is a simple demonstration program that calls EDWav.m on the Matlab platform

clear 

% load a Matlab formatted test seismic waveform data
load('Demo_data.mat');

% evaluate this waveform
obj=EDWav.evaluate(wvfm,Fs);

% plot the results
obj.result;

% You can add many restrictions to the evaluation, for example, the
% frequency band of the evaluation is set to [0 220]. 
band=[0 220]; % Since 'wvfm' at greater than 220Hz has been filtered, the effective frequency band is [0 220]
obj=EDWav.evaluate(wvfm,Fs,'FB',band);         % perform evaluation
obj.result;

% You can also adjust the frequency resolution of the results by setting
% the number of sampling points in the frequency domain.
FR=220; % frequency resolution
obj=EDWav.evaluate(wvfm,Fs,'FB',band,'FR',FR); % perform evaluation
obj.result;

% For more parameter settings, see the description in 'EDWav.evaluate'.

%%
% Here, we demonstrate the use of the 'EDWav.evaluate' function. While the 'EDWav.evaluate_sliding' function is not presented.
% The usage of the 'EDWav.evaluate_sliding' function can refer to the instructions in /Python Version/Demo_for_selection.ipynb'.

