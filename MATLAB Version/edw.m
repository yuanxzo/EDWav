% MIT License
% 
% Copyright (c) 2022
%     Bo YANG     (seism.yang@foxmail.com),
%     Haoran MENG (menghr@sustech.edu.cn) &
%     Ning GU     (guning@sustech.edu.cn)
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%
% =========================================================================
% EDWav Method
% This is a MATLAB package of the Advanced Evaluation Method for Diffuse
% Waveform (EDWav Method. Latest version https://github.com/yuanxzo/EDWav). 
%
%   This package contains the following static methods:
%     * obj=edw.edws;
%       obj=edw.read_sac;
%       obj=edw.read_mseed;
%       obj=edw.build_obj;
%
%   and normal methods:
%     * out=obj.evaluation;
%       out=obj.bandpass
%       out=obj.segment
%       out=obj.cut
%       out=obj.time
%       out=obj.freq
%       out=obj.fft
%       out=obj.find
%         H=obj.plot
%         H=obj.map
%     *   H=obj.result
%
%   Tips: MATLAB R2022a and later versions are recommended. The complete
%   operation of this package also requires the following additional
%   private methods, which users can request from the author of the package
%   by e-mail or download them from the links below.
%       load_sac.m        https://github.com/seismo-netizen/RN_detection/blob/master/load_sac.m
%       ReadMSEEDFast.m   https://www.mathworks.com/matlabcentral/fileexchange/46532-readmseedfast
%       progressBar.m     https://github.com/yuanxzo/Matlab-progressBar
%   
%    
%
%   Date:   2022.12.07 (Version 1.0.0)
%           2023.03.19 (Version 1.0.1)
%           2023.06.19 (Version 1.0.2)
%   References: 
%       [1] Bo Yang, Haoran Meng, Ning Gu, Xin Liu, Xiaofei Chen, and
%           Yehuda Ben-Zion. 2023. A Frequency Domain Methodology for
%           Quantitative Evaluation of Diffuse Wavefield with Applications
%           to Seismic Imaging. In preparation.
% =========================================================================

classdef edw
%% Each object of 'edw' has the following properties
properties
    data     % waveform data, the vector of npts*1;
    station  % some information about the station that records the 'data';
    A        % the condition A, which is a vector;
    B        % the condition B, which is a matrix;
    C        % the condition C, which is a matrix;
    F        % the frequency coordinate of A, B and C, which is a vector;
    proxy    % the proxies of A, B and C, which is a 3*1 vector;
    info     % some information about A, B and C
    isdiff   % The judgment of the evaluation result. If isdiff=1, it means that the evaluated data meets the requirements of diffuse field, otherwise not.
end

%% Static methods, which can be called for externally by edw.methodname
methods (Static)

    function obj=edws(data,varargin)
        % This function is the primary program provided by this package. It
        % can evaluate whether the waveforms contained in 'data' (a vector
        % of npts * 1) meet the conditions of a diffuse field. There are
        % many optional parameters for users to adjust according to
        % different problems to be handled. Enter name-value pair to
        % control a specified optional parameter, such as
        % obj=edw.edws(data,'NT',2), which means that two tapers are used
        % instead of the default one when evaluating.

        obj=edw;if isempty(data);return;end
        % Extract input parameters
        p = inputParser;% The following are optional parameters. When the input is not specified, each parameter uses its default value.
        addOptional(p,'Fs',1)    % sampling rate
        addOptional(p,'FB',[])   % frequency band of the data to be evaluated
        addOptional(p,'FR',100)  % frequency resolution of conditions A, B and C
        addOptional(p,'NT',1)    % number of tapers
        addOptional(p,'LT',[])   % length of taper
        addOptional(p,'TO',0.05) % tolerance of proxies of conditions A, B and C
        addOptional(p,'SE',[1,length(data(:,1))])  % indices of the starting and ending points of the data to be evaluated
        addOptional(p,'SF',0.05) % scale factor of the scale dependent RMS
        addOptional(p,'CP','no') % whether to perform an auxiliary calculation for finding the corner point of scale factor?
        addOptional(p,'PB','no') % show calculation progress? Tip: the display of work progress will increase the calculation time to a certain extent
        addOptional(p,'PL','no') % are parallel loops used? Tip: when the data to be evaluated is short, using parallelism does not necessarily reduce the calculation time
        p.parse(varargin{:});
        
                                        % Limitations (as follows)
        Fs              = p.Results.Fs;
        freq_band       = p.Results.FB;
        num_of_taper    = p.Results.NT;
        len_of_taper    = p.Results.LT;
        tol             = p.Results.TO; % 0<=tol<=1
        freq_res        = p.Results.FR; % 0<freq_res<?
        start_end_point = p.Results.SE;
        scale_factor    = p.Results.SF; % 0<=scale_factor<=1
        CP              = p.Results.CP;
        PB              = p.Results.PB;
        PL              = p.Results.PL; % 0<=int8(PL)<=maxNumCompThreads or PL='yes'/'no'

        % Create an object for the input data
        obj=edw.build_obj(data,'Fs',Fs);
        
        % Call obj.evaluation function for the created object to work out
        % the evaluation result.
        obj=obj.evaluation('Fs',Fs,'FB',freq_band,'NT',num_of_taper,...
                           'LT',len_of_taper,'SF',scale_factor,'TO',tol,...
                           'FR',freq_res,'SE',start_end_point,'CP',CP,...
                           'PB',PB,'PL',PL);
    end
    

    function obj = read_sac(filename,varargin)
        % Use the command: edw.read_sac(filename) to read one or more SAC 
        % files. The output is the 'edw' object(s).
        % Example: filename='name1.SAC' or filename={'name1.SAC','name2.SAC'}, then
        %          obj=edw.read_sac(filename) 
        %       or obj=edw.read_sac(filename,'npts',10)
        %       or obj=edw.read_sac(filename,'path','/mnt/d/folder/','npts',10)
        %
        %       or path='/mnt/d/folder/'; filename={dir([path,'*.SAC']).name};
        %          obj=edw.read_sac(filename,'path',path);

        % Extract input parameters
        p = inputParser;
        addOptional(p,'path',[])
        addOptional(p,'npts',[])
        addOptional(p,'PL','no')
        p.parse(varargin{:});
        
        path = p.Results.path;
        npts = p.Results.npts;        
        PL   = p.Results.PL;
        if strcmp(PL,'yes')
            parallel_threads=maxNumCompThreads;
        elseif isnumeric(PL)  %  0<=int8(PL)<=maxNumCompThreads
            PL=int8(PL);parallel_threads=PL;
            if PL<=0; parallel_threads=0;end
            if PL>maxNumCompThreads; parallel_threads=maxNumCompThreads;end
        else
            parallel_threads=0;
        end

        % Check the correctness of 'files' 
        if ischar(filename)
            N=1;filename={filename};
        elseif iscell(filename)
            N=length(filename);
        else
            error('There is an error in the input of files!')
        end

        % Read the sac file(s) and and create an object for each file
        obj(1:N)=edw;
        pbar = progressBar(N*2+1,'pname','READ SAC                  ');
        parfor (i=1:N, parallel_threads)
            obj(i)=edw;
            [hdr,data]=load_sac([path,filename{i}]);  % Quoted from: https://github.com/seismo-netizen/RN_detection/blob/master/load_sac.m
            
            obj(i).station.name=hdr.kstnm;
            obj(i).station.date=[num2str(hdr.nzyear),'J',num2str(hdr.nzjday)];
            obj(i).station.Fs=(1/hdr.delta);
            obj(i).station.dt=1/obj(i).station.Fs;
            obj(i).station.lat=hdr.stla;
            obj(i).station.lon=hdr.stlo;
            obj(i).station.npts=hdr.npts; 
            obj(i).data(1:obj(i).station.npts,1)=data;
            pbar.progress; %#ok<PFBNS>
        end
        
        % If 'npts' is input, which means the user requests that the length
        % of the data read from the sac file be equal to npts. If data is
        % shorter than npts, it will be zero-filled, while the part longer
        % than npts will be cut.
        pbar.pname='READ SAC => ALIGNMENT NPTS';pbar.progress;
        if ~isempty(npts)
            for i=1:N
                obj(i)=obj(i).cut([1,npts]);
                pbar.progress;
            end
        else
            npts=obj(1).station.npts;
            for i=1:N
                obj(i)=obj(i).cut([1,npts]);
                pbar.progress;
            end
        end
        pbar.stop;
    end
    

    function obj = read_mseed(filename,varargin)
        % Use the command: edw.read_mseed(filename) to read one or more
        % miniSEED files. The output is the 'edw' object(s).
        % Example:  filename='name1.mseed' or filename={'name1.mseed','name2.mseed'}.
        %           obj=edw.read_mseed(filename) 
        %        or obj=edw.read_mseed(filename,'npts',10)
        %        or obj=edw.read_mseed(filename,'path','/mnt/d/folder/','npts',10)
        %
        %        or path='/mnt/d/folder/'; filename={dir([path,'*.mseed']).name};
        %           obj=edw.read_mseed(filename,'path',path);

        % Extract input parameters
        p = inputParser;
        addOptional(p,'path',[])
        addOptional(p,'npts',[])
        addOptional(p,'lat',[])
        addOptional(p,'lon',[])
        addOptional(p,'PL','no')
        p.parse(varargin{:});
        
        path = p.Results.path;
        npts = p.Results.npts;
        lat  = p.Results.lat;
        lon  = p.Results.lon;
        PL   = p.Results.PL; 
        if strcmp(PL,'yes')
            parallel_threads=maxNumCompThreads;
        elseif isnumeric(PL)  %  0<=int8(PL)<=maxNumCompThreads
            PL=int8(PL);parallel_threads=PL;
            if PL<=0; parallel_threads=0;end
            if PL>maxNumCompThreads; parallel_threads=maxNumCompThreads;end
        else
            parallel_threads=0;
        end

        % Check the correctness of 'files' 
        if ischar(filename)
            N=1;filename={filename};
        elseif iscell(filename)
            N=length(filename);
        else
            error('There is an error in the input of files!')
        end
        
        if length(lat)~=length(lon)
            error('Longitude and latitude mismatch!')
        end
        if isempty(lon)
            lon=ones(N,1)*NaN;
        elseif ~isempty(lon) && length(lon) < N
            lon(length(lon)+1:N)=lon(length(lon));
        end
        if isempty(lat)
            lat=ones(N,1)*NaN;
        elseif ~isempty(lat) && length(lat) < N
            lat(length(lat)+1:N)=lat(length(lat));
        end
        
        % Read the miniSEED file(s) and and create an object for each file
        obj(1:N)=edw;
        pbar = progressBar(N*2+1,'pname','READ miniSEED                  ');
        parfor (i=1:N, parallel_threads)   % parfor may use more memory than 'for' loops. If out of memory, you can choose to set 'PR' to 'no'.
            obj(i)=edw;
            signalStruct=ReadMSEEDFast([path,filename{i}]);  % Quoted from: https://www.mathworks.com/matlabcentral/fileexchange/46532-readmseedfast
            dd=cat(1,signalStruct.data);
            
            obj(i).station.name=signalStruct(1).station;
            obj(i).station.date=[num2str(signalStruct(1).dateTime.year),'-',num2str(signalStruct(1).dateTime.month),'-',num2str(signalStruct(1).dateTime.day)];
            obj(i).station.Fs=signalStruct(1).sampleRate;
            obj(i).station.dt=1/obj(i).station.Fs;
            obj(i).station.lat=lat(i);
            obj(i).station.lon=lon(i);
            obj(i).station.npts=length(dd);
            obj(i).data(1:obj(i).station.npts,1)=dd(:);
            pbar.progress; %#ok<PFBNS>
        end
        
        % If 'npts' is input, which means the user requests that the length
        % of the data read from the miniSEED file be equal to npts. If data
        % is shorter than npts, it will be zero-filled, while the part
        % longer than npts will be cut.
        pbar.pname='READ miniSEED => ALIGNMENT NPTS';pbar.progress;
        if ~isempty(npts)
            for i=1:N
                obj(i)=obj(i).cut([1,npts]);
                pbar.progress;
            end
        else
            npts=obj(1).station.npts;
            for i=1:N
                obj(i)=obj(i).cut([1,npts]);
                pbar.progress;
            end
        end
        pbar.stop;
    end
    
    
    function obj=build_obj(data,varargin)  
        % Use the command: edw.build_obj(data) to build an 'obj'. 'data' is
        % a parameter that must be input. 
        % When 'Fs' is not input, it defaults to 1.
        % Example:  obj=edw.build_obj(data) 
        %        or obj=edw.build_obj(data,'Fs',500)
        %        or obj=edw.build_obj(data,'Fs',500,'name',{'R1'})

        if isempty(data);return;end
        [M,N]=size(data);

        % Extract input parameters
        Name(1:N,1)={'Unnamed'};
        Lat(1:N,1)={[]};
        Lon(1:N,1)={[]};
        p = inputParser;
        addOptional(p,'Fs',1)
        addOptional(p,'name',Name(1:N,1))
        addOptional(p,'lat',Lat)
        addOptional(p,'lon',Lon)
        p.parse(varargin{:});

        % Create an object for each column of data and write the parameters
        % to that object.
        obj(1:N)=edw;
        for i=1:N
            obj(i).data(1:M,1)  = data(:,i);
            if iscell(p.Results.name(i))
                obj(i).station.name = p.Results.name{i};
            else
                obj(i).station.name = p.Results.name(i);
            end
            obj(i).station.Fs   = p.Results.Fs;
            obj(i).station.dt   = 1/obj(i).station.Fs;
            obj(i).station.npts = M;
            if iscell(p.Results.lat(i))
                obj(i).station.lat  = p.Results.lat{i};
            else
                obj(i).station.lat  = p.Results.lat(i);
            end
            if iscell(p.Results.lon(i))
                obj(i).station.lon  = p.Results.lon{i};
            else
                obj(i).station.lon  = p.Results.lon(i);
            end
        end
    end
end


%% Normal methods, which can be used to operate on a created edw object.
methods
    
    function new_obj=evaluation(obj,varargin)
        % This function is called by edw.edws.
        % In fact, it can be used independently in the evaluation of the
        % data of an established object.
        % edw.edws is data oriented, but obj.evaluation is object oriented.
        % Therefore, edw.edws has the same effect as obj.evaluation. Users
        % can choose to use edw.edws or obj.evaluation to evaluate the data
        % they want to evaluate according to their needs.

        if isempty(obj.data);new_obj=obj;return;end
        % Extract input parameters
        p = inputParser;   % The optional parameters here are exactly the same as those in 'edw.edws'
        addOptional(p,'Fs',obj.station.Fs)
        addOptional(p,'FB',[])
        addOptional(p,'FR',100)
        addOptional(p,'NT',1)
        addOptional(p,'LT',[])
        addOptional(p,'TO',0.05)
        addOptional(p,'SE',[1,obj.station.npts])
        addOptional(p,'SF',0.05)
        addOptional(p,'CP','no')
        addOptional(p,'PB','no')
        addOptional(p,'PL','no')
        p.parse(varargin{:});

        Fs              = p.Results.Fs;
        freq_band       = p.Results.FB;
        num_of_taper    = p.Results.NT;
        len_of_taper    = p.Results.LT;
        scale_factor    = p.Results.SF;
        tol             = p.Results.TO;
        freq_res        = p.Results.FR;
        start_end_point = p.Results.SE;
        CP              = p.Results.CP;
        PB              = p.Results.PB;
        PL              = p.Results.PL;
        
        if isempty(freq_band)
            freq_band=[0,obj.station.Fs/2];
        end
        if strcmp(PL,'yes')
            parallel_threads=maxNumCompThreads;
        elseif isnumeric(PL)  %  0<=int8(PL)<=maxNumCompThreads
            PL=int8(PL);parallel_threads=PL;
            if PL<=0; parallel_threads=0;end
            if PL>maxNumCompThreads; parallel_threads=maxNumCompThreads;end
        else
            parallel_threads=0;
        end
        
        % Quality control of parameters
        min_len_of_taper=ceil(freq_res*Fs/(freq_band(2)-freq_band(1))); 
        if isempty(len_of_taper)
            len_of_taper=min_len_of_taper;
        else
            if len_of_taper<min_len_of_taper
                warning(['The length of the taper is too short, and the final evaluation result may be inaccurate. It is recommended to be greater than ', num2str(min_len_of_taper),'.']);
            end
        end
        
        num_of_win = floor((start_end_point(2)-start_end_point(1)+1)/len_of_taper); 
        if num_of_win<30
            warning(['The length of the data(:,1) is too short, and the final evaluation result may be inaccurate. It is recommended to be greater than ', num2str(30*len_of_taper),'.']);
        end
        [frq,fid,nfid]=obtain_freq(Fs,len_of_taper,freq_band);

        % Configure data
        wvfm=obj.data(start_end_point(1):start_end_point(2),1);
        idx = 1:len_of_taper:num_of_win*len_of_taper;
        wvfm_seg = zeros(len_of_taper,num_of_win);
        for i=1:num_of_win
            wvfm_seg(1:len_of_taper,i) = wvfm(idx(i):idx(i)+len_of_taper-1);
        end
	[tapers,weight] = sinusoidal_tapers(len_of_taper,num_of_taper);
        fft_all  = zeros(nfid,num_of_win,num_of_taper);
        for k=1:num_of_taper
            tmp = fft((wvfm_seg.*tapers(:,k)));
            fft_all(:,:,k) = tmp(fid(1)+1:fid(2)+1,:); 
        end
        
        % Evaluation with multitaper spectrum analysis
        tA = zeros(nfid,1);
        tB = zeros(nfid,nfid);
        tC = zeros(nfid,nfid);
        if strcmp(PB,'yes')
            pbar=progressBar(num_of_taper*num_of_win,'pname','EVALUATION');
        else
            pbar=progressBar(-1);
        end
        for k=1:num_of_taper
            wvfm_fft = fft_all(:,:,k); 
            E_power  = mean(abs(wvfm_fft).^2,2);
            [tmp_A,tmp_B,tmp_C,pbar] = calculation(wvfm_fft,E_power,nfid,num_of_win,pbar,parallel_threads);

            tA = tA+weight(k).*tmp_A;
            tB = tB+weight(k).*tmp_B;
            tC = tC+weight(k).*tmp_C;
        end
        pbar.stop;
        
        % Suppressing side lobe effect of sinusoidal tapers
        % (non final scheme)
        Cw=ones(nfid,nfid);CC=logspace(1,0,num_of_taper+1);
        for k=2:num_of_taper+1
            Cw(k:nfid+1:end)=CC(k-1);
        end
        Cw=Cw.*Cw';
        tC=tC./Cw;

        % Output
        new_obj=edw; new_obj.data=obj.data; new_obj.station=obj.station;
        
        [tproxy,new_obj.info.SF] = calculate_proxy(tA,tB,tC-eye(nfid),scale_factor,CP,PB); % Calculate the proxies of conditions A, B and C.
        if strcmp(CP,'yes')
            new_obj.proxy  = tproxy;
            new_obj.isdiff = [];
        else
            [~,sfin]=min(abs(scale_factor-new_obj.info.SF));
            new_obj.proxy  = tproxy(:,sfin);
            new_obj.isdiff = new_obj.proxy(1,1)<=tol && new_obj.proxy(2,1)<=tol && new_obj.proxy(3,1)<=tol;
        end
        
        temp = obj.bandpass(frq([fid(1),fid(2)]),'SE',start_end_point);
        new_obj.A = tA;
        new_obj.B = tB;
        new_obj.C = tC;
        new_obj.F = frq(fid(1):fid(2));
        new_obj.info.len_of_taper   = len_of_taper;
        new_obj.info.num_of_win     = num_of_win;
        new_obj.info.tapers         = tapers;
        new_obj.info.evaluated_data = temp.data;
        new_obj.info.band           = freq_band;
        new_obj.info.tol            = tol;
    end
    
    
    function new_obj=bandpass(obj,band,varargin)    
        % Perform band-pass filtering on an obj.data
        % Example:  band=[Fs/4, Fs/2]; id=3;
        %           new_obj=obj(id).bandpass(band);
        % Note:	if id=1 and length(obj)=1, id can be omitted, 
        %       i.e., new_obj=obj.bandpass(band);

        p = inputParser;
        addOptional(p,'SE',[1,obj.station.npts])
        p.parse(varargin{:});
        
        se = p.Results.SE;
        
        new_obj=obj;
        new_obj.station.npts=se(2)-se(1)+1;
        new_obj.data=obj.data(se(1):se(2),1);
        new_obj.info.band=band;
        Fs=obj.station.Fs;
        if isempty(band)
            warning('No input frequency band')
            return
        else
            if band(1)==0        % lowpass,  e.g., band=[0 Fs/2];
                [pB,pA]=butter(2,band(2)/(Fs/2),"low");
            elseif band(2)==Inf  % highpass, e.g., band=[Fs/2 Inf];
                [pB,pA]=butter(2,band(1)/(Fs/2),"high");
            else                 % bandpass, e.g., band=[Fs/4, Fs/2];
                [pB,pA]=butter(2,band/(Fs/2),'bandpass');
            end
            new_obj.data=filtfilt(pB,pA,obj.data(se(1):se(2),1));
        end
    end
    

    function new_obj=segment(obj,nseg,varargin)
        p = inputParser;
        addOptional(p,'SE',[1,obj(1).station.npts])
        p.parse(varargin{:});
        
        se  =p.Results.SE;
        if isempty(nseg)
            new_obj=obj;return;
        elseif nseg==1
            new_obj=obj.cut(se);return;
        else
            new_obj(nseg)=edw;
        end

        len=(se(2)-se(1)+1)/nseg;
        if rem(len,1)~=0
            error('Error input: NS')
        end
        odata=reshape(obj.data(se(1):se(2)),len,nseg);
        for i=1:nseg
            new_obj(i)=edw.build_obj(odata(:,i),'Fs',obj.station.Fs,'lat',obj.station.lat,'lon',obj.station.lon,'name',obj.station.name);
        end
    end


    function new_obj=cut(obj,varargin)
        % Perform data cutting on an obj.data with 'SE' 
        % Example: id=3; new_obj=obj(id).cut([1 100]);
        % Note:	if id=1 and length(obj)=1, id can be omitted, 
        %       i.e., new_obj=obj.cut([1 100]);

        p = inputParser;
        addOptional(p,'SE',[1,obj(1).station.npts])
        p.parse(varargin{:});
        se = p.Results.SE;
        
        new_obj=edw;
        new_obj.station=obj.station;
        
        tdata=obj.data;
        if new_obj.station.npts<se(2)
            tdata=[tdata;zeros(se(2)-obj.station.npts,1)];
        end
        new_obj.data=tdata(se(1):se(2),1);
        new_obj.station.npts=se(2)-se(1)+1;
    end
    

    function t=time(obj,varargin)
        % Get the time series of an 'obj' specified by 'id'.
        % Example:  id=3; t=obj(id).time;
        % Note:	if id=1 and length(obj)=1, id can be omitted, 
        %       i.e., t=obj.time;

        % Extract input parameters
        p = inputParser;
        addOptional(p,'SE',[1,obj.station.npts])
        p.parse(varargin{:});

        se = p.Results.SE;
        t(:,1)=(se(1)-1:se(2)-1).*obj.station.dt;
    end
    

    function f=freq(obj,varargin)
        % Get the frequency series of an 'obj' specified by 'id'.
        % Example:  id=3; f=obj(id).freq;
        % Note:	if id=1 and length(obj)=1, id can be omitted, 
        %       i.e., f=obj.freq;

        % Extract input parameters
        p = inputParser;
        addOptional(p,'SE',[1,obj.station.npts])
        p.parse(varargin{:});
        
        se = p.Results.SE;
        f(:,1)=(0:ceil(((se(2)-se(1))/2))).*(obj.station.Fs/(se(2)-se(1)+1));
    end  
    

    function spec=fft(obj,varargin)
        % Get the spectrum of fft of an 'obj.data'
        % Example:  id=3; spec=obj(id).fft; 
        %        or spec=obj(id).fft([100 10000]); 
        % Note:	if id=1 and length(obj)=1, id can be omitted, 
        %       i.e., spec=obj.fft;

        p = inputParser;
        addOptional(p,'SE',[1,obj.station.npts])
        p.parse(varargin{:});
        
        se = p.Results.SE;
        spec=fft(obj.data(se(1):se(2),1))/(se(2)-se(1)+1);
        spec=spec(1:ceil((se(2)-se(1)+2)/2),:);
        spec(2:end-1)=2*spec(2:end-1);
    end
   

    function id=find(obj,name)
        % Finds the ordinal number  mnvof the object with the specified station
        % name among all objects with its same variable name.
        % Example:  name='R0101'; id=obj.find(name); 

        if ischar(name)
            N=1;name={name};
        elseif iscell(name)
            N=length(name);
        else
            error('Please input the correct station name!')
        end
        
        id=zeros(N,1);
        for j=1:N
            nm=name{j};
            len=length(nm);
            for i=1:length(obj)
                if strcmp(obj(i).station.name(1:len),nm)
                    id(j)=i;
                    break
                end
            end
        end
    end
    
    
    function H=plot(obj,varargin)
        % Plot the normalized amplitude of data of object(s) in time and
        % frequency domain.
        % Example:  obj.plot or obj(2:10).plot([100,1000]);

        p = inputParser;
        addOptional(p,'SE',[])
        p.parse(varargin{:});
        
        se=p.Results.SE;
        N=length(obj);sname=cell(1,N);
        ntime=zeros(N,2);
        H=figure;
        subplot(1,2,1)
        for i=1:N
            if isempty(p.Results.SE)
                try
                    se=[1,obj(i).station.npts];
                catch 
                   continue
                end
            end
            time=obj(i).time([se(1),se(2)]);ntime(i,:)=[time(1),time(end)];
            plot(time,obj(i).data(se(1):se(2),1)./max(abs(obj(i).data(se(1):se(2),1)))+i*2);hold on
            sname{i}=obj(i).station.name;
        end
        xlabel('Time (s)');xlim([min(ntime(:,1)),max(ntime(:,2))])
        yticks(2:2:2*N);yticklabels(sname);ylim([1,2*N+1]);
        
        subplot(1,2,2)
        for i=1:N
            if isempty(p.Results.SE)
                se=[1,obj(i).station.npts];
            end
            temp=obj(i).fft([se(1),se(2)]);freq=obj(i).freq([se(1),se(2)]);
            semilogx(freq(2:end),abs(temp(2:end))./max(abs(temp(2:end)))+i);hold on
        end
        xlabel('Frequency (Hz)');
        yticks(1:N);yticklabels(sname);ylim([1,N+1]);
    end
    

    function H=map(obj)
        % Map the distribution of the stations recording data
        % Example:  id=1:10; obj(id).map;

        N=length(obj);
        lat=zeros(1,N);lon=zeros(1,N);name=cell(1,N);
        id=zeros(1,N);
        for i=1:N
            if isempty(obj(i).station) || isempty(obj(i).station.lat) || isempty(obj(i).station.lon)
                warning('These objects have no lat and lon attributes!')
                continue
            end
            id(i)=i;
            lat(i)=obj(i).station.lat;
            lon(i)=obj(i).station.lon;
            name{i}=['  ',obj(i).station.name];
        end
        if ~isempty(id(id>0))
            H=figure;
            scatter(lon(id),lat(id),'kv','filled','d');
            text(lon(id),lat(id),name(id),'fontsize',8)
            xlabel('longitude (Deg.)');ylabel('latitude (Deg.)')
        end
    end
    

    function H=result(obj)
        % When an obj.data has been evaluated, this function can be
        % used to draw images of conditions A, B and C.
        % Example:  id=3; obj(id).result;
        % Note:	if id=1 and length(obj)=1, id can be omitted, 
        %       i.e., obj.result;

        if isempty(obj.A) || isempty(obj.B) || isempty(obj.C)
            error('obj.data has not been evaluated!')
        end
        WW=23;HH=7;
        H=figure('Units','centimeters','Position',[15 7 WW HH]);
        subplot('Position',[1.5/WW,1.5/HH,4/WW,4/HH]);
        t=obj.time([1 length(obj.info.evaluated_data)]);
        plot(t,obj.info.evaluated_data./max(abs(obj.info.evaluated_data)),'k-');xlim([min(t) max(t)]);ylim([-1 1])
        xlabel('time (s)');axis square;ylabel('Amplitude (-)')
        title('(a) Normalized evaluated waveform_ ');
        % text(0.03*max(t),0.9,'(a)','fontsize',14)
        
        subplot('Position',[7/WW,1.5/HH,4/WW,4/HH]);
        plot(obj.F,obj.A,'k.-');xlim([min(obj.F) max(obj.F)]);ylim([0 1])
        xlabel('Frequency (Hz)');axis square
        title(['(b) Condition A     P_A = ',num2str(roundn(obj.proxy(1),-3))]);
        % text(0.03*(max(obj.F)-min(obj.F)),0.95,'(b)','fontsize',14)
        
        subplot('Position',[12.5/WW,1.5/HH,4/WW,4/HH]);
        imagesc(obj.F,obj.F,obj.B);clim([0 1])
        xlabel('Frequency (Hz)');ylabel('Frequency (Hz)');axis square
        title(['(c) Condition B     P_B = ',num2str(roundn(obj.proxy(2),-3))]);
        % text(0.85*(max(obj.F)-min(obj.F)),0.05*(max(obj.F)-min(obj.F)),'(c)','fontsize',14,'color','w')
        
        subplot('Position',[18/WW,1.5/HH,4/WW,4/HH]);
        imagesc(obj.F,obj.F,obj.C);clim([0 1])
        xlabel('Frequency (Hz)');ylabel('Frequency (Hz)');axis square
        title(['(d) Condition C     P_C = ',num2str(roundn(obj.proxy(3),-3))]);
        % text(0.85*(max(obj.F)-min(obj.F)),0.05*(max(obj.F)-min(obj.F)),'(d)','fontsize',14,'color','w')
    end
end
end

%% Private functions, which are the support functions for the 'edw' package

% Calculate conditions A, B and C
function [tmp_A,tmp_B,tmp_C,pbar] = calculation(wvfm_fft,E_power,nfrq,nwin,pbar,parallel_threads)
    Emn      = zeros(nfrq,nfrq);
    tmp_B    = zeros(nfrq,nfrq);
    Emn_str  = zeros(nfrq,nfrq);
    tmp_C    = zeros(nfrq,nfrq);
    
    parfor (i=1:nwin,parallel_threads)
        [emn,emn_str]=sub(wvfm_fft,i,nfrq);
        Emn     = Emn+emn./nwin;
        Emn_str = Emn_str+emn_str./nwin;
        pbar.progress; %#ok<PFBNS>
    end

    tmp_A = abs(mean(wvfm_fft,2)).^2./E_power;
    tmpE=E_power(:);
    for j=1:nfrq
        tmp_B(j,:) = abs(Emn(j,:)).^2./(E_power(j).*tmpE)';
        tmp_C(j,:) = abs(Emn_str(j,:)).^2./(E_power(j).*tmpE)';
    end 
end
function [emn,emn_str] = sub(wvfm_fft,i,nfrq)
    emn      = zeros(nfrq,nfrq);
    emn_str  = zeros(nfrq,nfrq);
    for j=1:nfrq
        emn(j,j:end) = wvfm_fft(j,i).*wvfm_fft(j:end,i);
        emn(:,j)     = emn(j,:);
        emn_str(j,j:end) = wvfm_fft(j,i).*conj(wvfm_fft(j:end,i));
        emn_str(:,j)     = conj(emn_str(j,:));
    end
end

% Calculate the scale dependent RMS of A, B and C
function [proxy,SF]=calculate_proxy(A,B,C,SF,CP,PB)
	[m,~]=size(B);
    if strcmp(CP,'yes') || strcmp(PB,'YES') || strcmp(CP,'YES')
		SS=1:m;SF=SS./m;
    else
        SS=ceil(m.*SF);
    end
    N=length(SS);
    if N>1; pbar = progressBar(N,'pname','SRMS'); end
    if N==1 && (strcmp(PB,'yes') || strcmp(PB,'YES')); pbar = progressBar(m,'pname','SRMS'); end
	proxy=zeros(3,N);
    for is=1:N
		S=SS(is);  WA=zeros(m,1);WB=zeros(m,m);WC=zeros(m,m);
        if S==0
            WA=A;WB=B;WC=C;
        elseif S==m
            WA=mean(A);WB=mean(B(:));WC=mean(C(:));
        else
            tid1=0;tid2=0;
            for i=1:m
                id1=max(1,i-S);id2=min(m,i+S);nid=id2-id1+1;
                if id1==tid1 && id2==tid2
                    WA(i)=WA(i-1);
                else
                    WA(i)=mean(A(id1:id2));
                end
                tid1=id1;tid2=id2;
                
                tjd1=0;tjd2=0;
                for j=i:m
                    jd1=max(1,j-S);jd2=min(m,j+S);njd=jd2-jd1+1;
                    if jd1==tjd1 && jd2==tjd2
                        WB(i,j)=WB(i,j-1);WC(i,j)=WC(i,j-1);
                    else
                        WB(i,j)=sum(sum(B(id1:id2,jd1:jd2)))/(njd*nid);
                        WC(i,j)=sum(sum(C(id1:id2,jd1:jd2)))/(njd*nid);
                    end
                    WB(j,i)=WB(i,j); WC(j,i)=WC(i,j);
                    tjd1=jd1;tjd2=jd2;
                end

                if N==1 && (strcmp(PB,'yes') || strcmp(PB,'YES'))
                    pbar.progress;
                end
            end
        end
		
		AS=A.*WA./mean(A);
		BS=B.*WB./mean(B(:));
		CS=C.*WC./mean(C(:));
		proxy(1,is)=sqrt(sum(AS(:).^2)/m);
		proxy(2,is)=sqrt(sum(BS(:).^2)/m^2);
		proxy(3,is)=sqrt(sum(CS(:).^2)/m^2);
        
        if N>1
            pbar.progress;
        end
    end
    if N>1 || strcmp(PB,'yes') || strcmp(PB,'YES')
        pbar.stop;
    end
        
	if strcmp(CP,'YES')
		figure;hold on
		plot(SF(1:end-1),diff(proxy(1,:)./max(abs(proxy(1,:)))),'r.-');
		plot(SF(1:end-1),diff(proxy(2,:)./max(abs(proxy(2,:)))),'g.-');
		plot(SF(1:end-1),diff(proxy(3,:)./max(abs(proxy(3,:)))),'b.-');
		xlabel('Scale factor');ylabel('First derivative of normalized SRMS')
		legend('Condition A','Condition B','Condition C');
		axis square;hold off
	end
end

% Sinusoidal tapers
function [tapers,lambda] = sinusoidal_tapers(len_of_win,num_of_taper)
    points=1:len_of_win;
    tapers=zeros(len_of_win,num_of_taper);
    for i=1:num_of_taper
        tapers(:,i)=sqrt(2/(len_of_win+1)).*sin((pi*i.*points)./(len_of_win+1));
    end
    lambda(1:num_of_taper)=1/num_of_taper;
end

% Obtain frequency information related to conditions A, B and C
function [frq,fid,nfid]=obtain_freq(Fs,len_of_taper,band)
    len_of_frq = ceil((len_of_taper-1)/2)-1;
    df         = Fs/len_of_taper;
    frq        = (1:1:len_of_frq).*df;  % Zero frequency is not considered

    fid=[1,len_of_frq];
    for i=1:len_of_frq-1
        if band(1)>frq(i) && band(1)<frq(i+1)
            fid(1)=i+1;
        end
        if band(1)==frq(i)
            fid(1)=i;
        end
        if band(2)>frq(i) && band(2)<frq(i+1)
            fid(2)=i;
        end
        if band(2)==frq(i+1)
            fid(2)=i+1;
        end
    end
    nfid=fid(2)-fid(1)+1;
end
