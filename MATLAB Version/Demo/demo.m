%% A demonstration of diffuse waveform evaluation by 'edw' package
% The 'edw' is a MATLAB package of the advanced evaluation method for 
% diffuse waveform (EDWav Method). It provides two methods (functions) to
% evaluate seismic waveform.

% Please run section by section.
%% First, the data oriented method: 'edw.edws' function.
% When you only have time series data instead of SAC or miniSEED files, you
% can use 'edw.edws' function to evaluate this data.
clear
load('demo_data.mat')
obj=edw.edws(wvfm,'Fs',Fs);                   % perform evaluation
obj.result;

% You can add many restrictions to the evaluation, for example, the
% frequency band of the evaluation is set to [0 220]. 
band=[0 220]; % Since 'wvfm' at greater than 220Hz has been filtered, the effective frequency band is [0 220]
obj=edw.edws(wvfm,'Fs',Fs,'FB',band);         % perform evaluation
obj.result;

% You can also adjust the frequency resolution of the results by setting
% the number of sampling points in the frequency domain.
FR=220; % frequency resolution
obj=edw.edws(wvfm,'Fs',Fs,'FB',band,'FR',FR); % perform evaluation
obj.result;

% For more parameter settings, see the description in 'edw.edws'.
%% The second method is object oriented: 'obj.evaluation' function.
% The precondition for using this method is that you already have an 'edw'
% object. There are two ways to create an edw object. The first is to use
% the function edw.build_obj, and the second is to use the edw.read_sac or
% edw.read_mseed functions.

% The first way to create an edw object is by edw.build_obj
clear
load('demo_data.mat')
obj1=edw.build_obj(wvfm,'Fs',Fs); 
obj1=obj1.evaluation;              % perform evaluation
obj1.result;

% The second way to create an edw object is by edw.read_sac or edw.read_mseed
obj2=edw.read_mseed({'ANMO.mseed','COLA.mseed'},'npts',600*40,'lat',[34.946 64.874],'lon',[-106.457 -147.862]); % Here, we additionally add a restriction on the length of the read data and assign each data its station position
for i=1:length(obj2)  % since two files are read in, you can use the for loop to evaluate all 'obj'
    obj2(i)=obj2(i).evaluation;   % perform evaluation
    obj2(i).result;
end

%% A brief introduction to some other functions of 'edw' package
obj2.plot;  % view the waveform in time and frequnecy domian
obj2.map;   % view the stations distribution

t    = obj1.time;  % calculate time, frequency series, and spectrum of an 'obj'
f    = obj1.freq;
spec = obj1.fft;

band=[1 10];    % band pass filtering before evaluation
obj2_bp=obj2(1).bandpass(band);
obj2_bp=obj2_bp.evaluation('FB',band);
obj2_bp.result;

obj2_seg=obj2(1).segment(10,'SE',[1 600*40]); % segment the selected data into multiple objects

obj2_cut=obj2(1).cut('SE',[1 600*40]);  % keep selected data

%% How to select the appropriate scale factor of the scale dependent RMS?
% First, switch on the 'CP' and perform evaluation, and an image will appear as a result
temp=obj2_bp.evaluation('FB',band,'CP','YES'); 

% Then, judge the general corner point of the three curves with the naked
% eye, and the abscissa of this point is the value of 'SF' that we think is
% appropriate. 
% As you can see from the figure, the value corresponding to this obj is
% roughly 0.06. In our experience, the value of 'SF' is basically in the
% range of [0.03,0.07], so the default value of 'SF' we set in the program
% is 0.05. Generally, this default value is feasible, and users do not
% have to perform these steps every time.

% Finally, modify the value of 'SF', turn off 'CP', and re-evaluate.
SF=0.07;
new=obj2_bp.evaluation('FB',band,'SF',SF); 
new.result;
% Whether it is from the images of conditions A, B and C or from the
% calculated proxies of them (located on the title of the figure), it can
% be found that 'obj2_bp' satisfies the diffuse field characteristics (We
% declare that the object with PA, PB and PC all less than 0.05 conform to
% the diffuse field characteristics. 0.05 is the default value. Different
% values can be set for different data.), so the value of 'new.isdiff'
% is 1.
disp(new.isdiff)
%% At the end of demonstration
% More optional parameters are not introduced one by one, their meanings
% can be found in the description of optional parameters of 'edw.edws'. 
%
% Please look forward to our paper for the specific applications and
% physical implications of this method:
%   [1] Bo Yang, Haoran Meng, Ning Gu, Xin Liu, Xiaofei Chen, and Yehuda
%       Ben-Zion. 2023. A Frequency Domain Methodology for Quantitative
%       Evaluation of Diffuse Wavefield with Applications to Seismic
%       Imaging. Preprint.
