classdef EDWav
%EDWav Evaluation of the Diffuseness of the Wavefield Method
%   Main functions: 
%
%   Note: 
%
% Latest date: 2024-03-22
%
% Copyright (C) 2019 <Bo Yang. Email:seism.yang@foxmail.com>

properties
    A=0         % 描述弥散程度结果的条件A
    B=0         % 描述弥散程度结果的条件B
    C=0         % 描述弥散程度结果的条件C
    proxy=[]    % 描述弥散程度结果的指标
    freqs=[]    % 弥散程度结果的依赖频率
    stats=[]    % 与评估结果相关的一些信息
end

methods (Static = true)
    function obj=evaluate(wave,FS,options)
        % 函数功能介绍: 评估输入地震波形段的弥散程度
        arguments
            wave (:,1) {mustBeNumeric}  % 单道地震波形数据, size(wave) is [npts,1]
            FS  (1,1) {mustBeNumeric}=1 % 波形数据采样率 sampling rate
            options.FB (1,2) {mustBeNumeric} =[0 FS/2] % 要评估的频率范围, 例如: [1, 100], 单位Hz;
            options.FR (1,1) {mustBeNumeric} =100 % 结果的最小频率分辨率, 即FB被等分的最少个数 frequency resolution of conditions A, B and C
            options.NT (1,1) {mustBeNumeric} =1  % Multitaper的个数 umber of tapers
            options.SF (1,1) {mustBeNumeric} =0.05 % sRMS的一个参数.scale factor of the scale dependent RMS
        end
        FB = options.FB;
        FR = options.FR;
        NT = options.NT;
        SF = options.SF;

        obj=EDWav;
        NT =ceil(NT);
        FR =ceil(FR);

        % 波形采样点数
        npts=length(wave(:,1)); 

        % 要解析的频率范围
        if isempty(FB)
            FB=[0, FS/2];
        end
        
        % 根据输入的FR,FB计算窗口的最小长度
        len_of_win=ceil(FR*FS/(FB(2)-FB(1))) ;
        if npts < len_of_win
            error("波形长度过小, 结束评估!");
        end
        
        % 总波形被该长度的窗分割的个数
        num_of_win=floor(npts/len_of_win);    
        if num_of_win < 30
            disp("   警告: 根据输入的设置, 最终结果可能不准确. 建议输入的波形长度 npts>" + num2str(ceil(30*len_of_win)))
        end

        % 最终窗口的长度
        len_of_win=floor(len_of_win+(npts-num_of_win*len_of_win)/num_of_win);

        % 该窗口能解析的频率
        freqs=0 : FS/len_of_win : FS/2+FS/len_of_win;
        if freqs(end)>FS/2
            freqs(end)=[];
        end

        % 需要解析的频率
        fin = freqs>=FB(1) & freqs<=FB(2);
        obj.freqs=freqs(fin);
        nfin=length(obj.freqs);

        % 加taper
        [tapers, weight] = sinusoidal_tapers(len_of_win,NT);

        % 分割波形
        wave_seg=zeros(len_of_win,num_of_win);
        for i = 1:num_of_win
            wave_seg(:,i)=wave((i-1)*len_of_win+1 : i*len_of_win);
        end

        % 滤波到指定频带
        % if FB(1)~=0 && (FB(2)~=Inf || FB(2)~=FS/2)
        %     if FB(1)==0        % lowpass,  e.g., band=[0 FS/2];
        %         [pB,pA]=butter(2,FB(2)/(FS/2),"low");
        %     elseif FB(2)==Inf || FB(2)==FS/2  % highpass, e.g., band=[FS/4 Inf];
        %         [pB,pA]=butter(2,FB(1)/(FS/2),"high");
        %     else               % bandpass, e.g., band=[FS/4, FS/2];
        %         [pB,pA]=butter(2,FB   /(FS/2),'bandpass');
        %     end
        %     wave_seg=filtfilt(pB,pA,wave_seg);
        % end
        

        % 开始计算弥散场的三个条件A, B, C
        for i = 1:NT
            wave_fft = fft(wave_seg.*tapers(:,i), [], 1)/len_of_win;
            wave_fft = wave_fft(1:ceil(len_of_win/2+1),:);
            wave_fft(2:end-1,:) = 2*wave_fft(2:end-1,:);
            wave_fft = wave_fft(fin,:);

            E_power = mean(abs(wave_fft).^2, 2);
            EEpower = E_power .* E_power';

            A = abs(mean(wave_fft, 2)) .^2 ./ E_power;
            tempB = 0; tempC = 0;
            for j = 1 : num_of_win
                tempB = tempB + wave_fft(:,j) .* wave_fft(:,j).';
                tempC = tempC + wave_fft(:,j) .* conj(wave_fft(:,j).');
            end
            B = abs(tempB/num_of_win) .^2 ./ EEpower;
            C = abs(tempC/num_of_win) .^2 ./ EEpower;

            obj.A = obj.A + weight(i).*A;
            obj.B = obj.B + weight(i).*B;
            obj.C = obj.C + weight(i).*C;
        end

        % 压制条件C的旁瓣
        Cw=ones(nfin,nfin);CC=logspace(1,0,NT+1);
        for k=2:NT+1
            Cw(k:nfin+1:end)=CC(k-1);
        end
        Cw=Cw.*Cw';
        obj.C=obj.C./Cw;


        % 计算三个条件的弥散度大小
        if SF>=0 && SF<=1
            obj = obj.sRMS(SF);
        else
            obj.proxy=[NaN,NaN,NaN];
        end

        % 保留一些与结果相关的信息
        obj.stats.len_of_win = len_of_win;
        obj.stats.num_of_win = num_of_win;
        obj.stats.num_of_taper = NT;
        obj.stats.SF_of_sRMS = SF;
        obj.stats.FB = FB;

    end

    function [proxy,stats]=evaluate_sliding(wave,FS,options)
        % 函数功能介绍: 以滑动窗评估连续地震波形（通常为小时级别或更久的波形）的弥散程度
        % 输出参数介绍: proxy: 一个弥散度随时间变化的值; stats: 与计算相关的一些信息量
        arguments
            wave (:,1) {mustBeNumeric}  % 单道地震波形数据, size(wave) is [npts,1]
            FS  (1,1) {mustBeNumeric}=1 % 波形数据采样率 sampling rate
            options.FB (1,2) {mustBeNumeric} =[] % 要评估的频率范围, 例如: [1, 100], 单位Hz;
            options.FR (1,1) {mustBeNumeric} =100 % 结果的最小频率分辨率, 即FB被等分的最少个数 frequency resolution of conditions A, B and C
            options.NT (1,1) {mustBeNumeric} =1  % Multitaper的个数 umber of tapers
            options.SF (1,1) {mustBeNumeric} =0.05 % sRMS的一个参数.scale factor of the scale dependent RMS
            options.NW (1,1) {mustBeNumeric} =30 % 滑动窗所包含的小窗的最少个数, 建议该值最小为30, 通常取值在[30,60]之间.
            options.DW (1,1) {mustBeNumeric} =30  % 滑动窗的间隔小窗数
            options.PN ='' % progressBar's name 
        end
        FB = options.FB;
        FR = options.FR;
        NT = options.NT;
        SF = options.SF;
        NW = options.NW;
        DW = options.DW;

        obj=EDWav;
        NT =ceil(NT);
        FR =ceil(FR);
        NW =ceil(NW);
        DW =ceil(DW);

        % 波形采样点数  
        npts=length(wave(:,1)); 

        % 要解析的频率范围
        if length(FB) ~= 2
            FB=[0, FS/2];
        end
        
        % 根据输入的FR,FB计算小窗口的最小长度
        len_of_win=ceil(FR*FS/(FB(2)-FB(1))) ;
        if npts < len_of_win * NW
            error("滑动窗的长度已大于波形长度, 结束评估!");
        end
                 
        % 滑动窗内的小窗能解析的频率
        freqs=0 : FS/len_of_win : FS/2+FS/len_of_win;
        if freqs(end)>FS/2
            freqs(end)=[];
        end

        % 需要解析的频率
        fin= freqs>=FB(1) & freqs<=FB(2) ;
        obj.freqs=freqs(fin);
        nfin=length(obj.freqs);

        % 该连续波形允许的小窗的个数
        num_of_win=floor(npts/len_of_win);

        % 该连续波形允许的滑动窗的个数
        sindex=1:DW:num_of_win-NW+1;
        num_of_sld=length(sindex);
        
        % 最终小窗的长度
        len_of_win=floor(len_of_win+(npts-num_of_win*len_of_win)/num_of_win);

        % 加taper
        [tapers, weight] = sinusoidal_tapers(len_of_win,NT);

        % 分割波形
        wave_seg=zeros(len_of_win,num_of_win);
        for i = 1:num_of_win
            wave_seg(:,i)=wave((i-1)*len_of_win+1 : i*len_of_win);
        end
        
        % % 滤波到指定频带
        % if FB(1)~=0 && (FB(2)~=Inf || FB(2)~=FS/2)
        %     if FB(1)==0        % lowpass,  e.g., band=[0 FS/2];
        %         [pB,pA]=butter(3,FB(2)/(FS/2),"low");
        %     elseif FB(2)==Inf || FB(2)==FS/2  % highpass, e.g., band=[FS/4 Inf];
        %         [pB,pA]=butter(3,FB(1)/(FS/2),"high");
        %     else               % bandpass, e.g., band=[FS/4, FS/2];
        %         [pB,pA]=butter(3,FB   /(FS/2),'bandpass');
        %     end
        %     wave_seg=filtfilt(pB,pA,wave_seg);
        % end

        % 开始滑动计算弥散场的三个条件A, B, C
        wave_fft=zeros(length(obj.freqs),num_of_win,NT);
        for i = 1 : NT
            wave_tmp = fft(wave_seg.*tapers(:,i), [], 1)/len_of_win;
            wave_tmp = wave_tmp(1:ceil(len_of_win/2+1),:);
            wave_tmp(2:end-1,:) = 2*wave_tmp(2:end-1,:);
            wave_fft(:,:,i) = wave_tmp(fin,:);
        end

        proxy=zeros(num_of_win*len_of_win,1);
        stats.stack=zeros(num_of_win*len_of_win,1);

        % 开始滑动评估
        if strcmp(options.PN,'')
            NN=-1;
        else
            NN=num_of_sld;
        end

        % 分母
        E_power = cell(NT,1);
        EEpower = cell(NT,1);
        for i = 1:NT
            E_power{i} = median(abs(wave_fft(:,:,i)).^2, 2);
            EEpower{i} = E_power{i}.*E_power{i}';
        end

        pp=progressBar(NN,'pname',options.PN);
        for k = 1 : num_of_sld
            index=[sindex(k) sindex(k)+NW-1];
            indexs = (index(1)-1)*len_of_win+1 : index(2)*len_of_win;

            % 计算该滑动窗的弥散性条件
            obj.A=0; obj.B=0; obj.C=0;

            for i = 1:NT
                wave_use= wave_fft(:,index(1):index(2),i);
    
                A = abs(mean(wave_use, 2)) .^2 ./ E_power{i};
                tempB = 0; tempC = 0;
                for j = 1 : NW
                    tempB = tempB + wave_use(:,j) .* wave_use(:,j).';
                    tempC = tempC + wave_use(:,j) .* conj(wave_use(:,j).');
                end
                B = abs(tempB/NW) .^2 ./ EEpower{i};
                C = abs(tempC/NW) .^2 ./ EEpower{i};
    
                obj.A = obj.A + weight(i).*A;
                obj.B = obj.B + weight(i).*B;
                obj.C = obj.C + weight(i).*C;
            end

            % 判断结果中是否有NaN, 有的话即将该段结果置为1
            if any(isnan(obj.A)) || any(isnan(obj.B(:))) || any(isnan(obj.C(:)))
                proxy(indexs, 1)=1;
                stats.stack(indexs, 1)  =stats.stack(indexs, 1)+0;
            else
                % 压制条件C的旁瓣
                Cw=ones(nfin,nfin);CC=logspace(1,0,NT+1);
                for kk=2:NT+1
                    Cw(kk:nfin+1:end)=CC(kk-1);
                end
                Cw=Cw.*Cw';
                obj.C=obj.C./Cw;

                % 计算三个条件的弥散度大小
                obj = obj.sRMS(SF);

                % 记录结果
                proxy(indexs, 1)=proxy(indexs, 1)+mean(obj.proxy);  % 忽略obj.C
                stats.stack(indexs, 1)=stats.stack(indexs, 1)+1;
            end

            pp.progress;
        end
        pp.stop;            

        % 保留一些与结果相关的信息
        proxy=proxy./stats.stack;
        
        stats.freqs = obj.freqs;
        stats.len_of_win = len_of_win;
        stats.num_of_win = num_of_win;
        stats.num_of_sld = num_of_sld;
        stats.num_of_taper = NT;
        stats.SF_of_sRMS = SF;
        stats.NW = NW;
        stats.DW = DW;
        stats.FB = FB;

    end
end

methods
    function obj=sRMS(obj,SF)
        % 计算X的sRMS函数值
        N = floor(length(obj.A(:)));
        s = ceil(N*SF);

        obj_C=obj.C-eye(length(obj.A));

        meanA=mean(obj.A);
        meanB=mean(obj.B,"all");
        meanC=mean(obj_C,"all");
    
        wA=zeros(N,1);
        wB=zeros(N,N);
        wC=zeros(N,N);
    
        if s==N
            wA=1;
            wB=1;
            wC=1;
        else
            for i = 1 : N
                indexi = i-s : 1 : i+s;
                fin= indexi>0 & indexi<=N ;
                indexi = indexi(fin);

                wA(i)=mean(obj.A(indexi));

                for j = i : N
                    indexj = j-s : 1 : j+s;
                    fin= indexj>0 & indexj<=N ;
                    indexj = indexj(fin);

                    wB(i,j)=mean(obj.B( indexi, indexj ),"all");
                    wC(i,j)=mean(obj_C( indexi, indexj ),"all");
                    wB(j,i)=wB(i,j);
                    wC(j,i)=wC(i,j);
                end
            end

            wA=wA./meanA;
            wB=wB./meanB;
            wC=wC./meanC;
        end

        obj.proxy=zeros(1,3);
        obj.proxy(1)=sqrt(mean((wA .* obj.A) .^2,"all"));
        obj.proxy(2)=sqrt(mean((wB .* obj.B) .^2,"all"));
        obj.proxy(3)=sqrt(mean((wC .* obj_C) .^2,"all"));
        % Considering the similarity between conditions B and C, it is possible to make the proxy of condition C 
        % equal to that of condition B to avoid the influence of diagonal sidelobes on the proxy of condition C
        % (Note that this is only an empirical operation. You can choose to comment on the following line of code to reject this operation.)
        obj.proxy(3)=obj.proxy(2);
    end

    function H=result(obj)
        if isempty(obj.A) || isempty(obj.B) || isempty(obj.C)
            error('obj.data has not been evaluated!')
        end

        WW=17;HH=7;
        H=figure('Units','centimeters','Position',[15 10 WW HH]);

        subplot('Position',[1/WW,1.5/HH,4/WW,4/HH]);
        plot(obj.freqs,obj.A,'k.-');xlim([min(obj.freqs) max(obj.freqs)]);ylim([-0.1 1])
        xlabel('Freq. (Hz)');axis square
        title(['Cond. A     P_A = ',num2str(roundn(obj.proxy(1),-3))]);
        
        subplot('Position',[6.5/WW,1.5/HH,4/WW,4/HH]);
        imagesc(obj.freqs,obj.freqs,obj.B);clim([0 1])
        xlabel('Freq. (Hz)');ylabel('Freq. (Hz)');axis square
        title(['Cond. B     P_B = ',num2str(roundn(obj.proxy(2),-3))]);
        
        subplot('Position',[12/WW,1.5/HH,4/WW,4/HH]);
        imagesc(obj.freqs,obj.freqs,obj.C);clim([0 1])
        xlabel('Freq. (Hz)');ylabel('Freq. (Hz)');axis square
        title(['Cond. C     P_C = ',num2str(roundn(obj.proxy(3),-3))]);
    end

end

end


% Sinusoidal tapers
function [tapers,lambda] = sinusoidal_tapers(len_of_win,num_of_taper)
    points=1:len_of_win;
    tapers=zeros(len_of_win,num_of_taper);
    for i=1:num_of_taper
        tapers(:,i)=sqrt(2/(len_of_win+1)).*sin((pi*i.*points)./(len_of_win+1));
    end
    % lambda(1:num_of_taper)=1/num_of_taper;
    lambda(1:num_of_taper)=(num_of_taper:-1:1)./sum(num_of_taper:-1:1);
end
