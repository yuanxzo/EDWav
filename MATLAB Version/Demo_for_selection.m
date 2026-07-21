clear 
% read test wavefrom data by obspy
files=dir('../Python Version/Demo_data/*.new');
for i=1:length(files)  % There are 5 "*.sac.new" files in the folder "Demo_data".
    [sachdr,temp]=load_sac(['../Python Version/Demo_data/',files(i).name]); 
    Data(1:sachdr.npts,i)=temp;
end

% bandpass filtering
FS=1/sachdr.delta;
[pB,pA]=butter(4,[0.03 0.5]/(FS/2),'bandpass');
for i=1:length(Data(1,:))
    Data(:,i)=filtfilt(pB,pA,Data(:,i));
end

% Specify the waveform to perform the selection operation
waveform=Data(1728001:3456000,:); % datetime(2011,7,13,0,0,0) to datetime(2011,7,14,0,0,0);
[npts,nsta]=size(waveform);

% set the parameters for sliding evaluation
% -Sampling rate: FS.
% -Frequency band range for analysis: FB.
% -Frequency resolution: FR.
% -The length of the sliding window: NW. NW represents the length of the sliding window is NW times the length of the small window. The length of the small window is defined by "len_of_win" in the "evaluate" function and is regulated by "FR" and "FB". Generally, the length of a sliding window can be consistent with the length of subsequent cross-correlation operations.
% -The movement interval of the sliding window: DW. DW represents the movement interval of the sliding window is DW times the length of the small window.
FB=[0.03,0.5];
FR=47;
NW=6;
DW=6;  %DW=NW means no overlap between two adjacent sliding windows


% begin to evaluate the waveform data with sliding window
P=zeros(npts,nsta);  % P is to store the diffuseness proxy value of the waveform data. 1728000 is the length of the waveform data. 
tic
for i=1:nsta
    [proxy, stats] = EDWav.evaluate_sliding(waveform(:,i), FS, 'FB',FB, 'FR',FR, 'NW',NW, 'DW',DW);
    P(:,i)=proxy;  % store the diffuseness proxy
end
toc
disp(['The length of the small window is ',num2str(stats.len_of_win/FS),' s.'])
disp(['The length of the sliding window is ',num2str(NW*stats.len_of_win/FS),' s.'])


% Now, set a threshold to select the diffuse waves
% We can use the MAD (median absolute deviation) method to determine the threshold. The threshold is calculated by: threshold=3*MAD+median
% A threshold can be calculated for each waveform data or a uniform threshold can be calculated for all waveform data.
% Here, we calculate a threshold for each waveform data.
threshold=zeros(nsta,1);
P_selected=zeros(npts,nsta);
waveform_selected=zeros(npts,nsta);
for i=1:nsta
    % We can operate on the P value or the logarithmic P value, which depends on your data
    threshold(i)=3*median(abs(log(P(:,i))-median(log(P(:,i)),1,"omitmissing")),1,"omitmissing")+median(log(P(:,i)),1,"omitmissing");

    % select the diffuse waves, which are the waveform data with P value less than the threshold
    P_selected(:,i)=log(P(:,i)).*(log(P(:,i))<=threshold(i));
    P_selected(P_selected(:,i)==0,i)=NaN;
    waveform_selected(:,i)=waveform(:,i).*sign(abs(P_selected(:,i)));
end

% plot the selected waveform data
step = 20;
timeIdx = 1:step:npts; 
t = datetime(2011,7,13,0,0,0) + seconds((timeIdx-1)/FS);

figure('Position', [100, 100, 1200, 800]);
for i = 1:size(waveform, 2)   
    subplot(5, 1, i);
    
    % Left axis: Original waveform and selected waveform
    colororder(gca,{'k','k'})
    yyaxis left;
    plot(t, waveform(timeIdx, i),'-', 'Color', [0.75, 0.75, 0.75],'LineWidth',1.5); hold on;
    plot(t, waveform_selected(timeIdx, i),'-','Color',[0.30,0.75,0.93],'LineWidth',1.5);
    ylabel(files(i).name(1:4));
    yticklabels([]);
    
    % Right axis: P-value (original and selected) and threshold line
    yyaxis right;
    plot(t, log(P(timeIdx, i)), 'k-','LineWidth',1.5);
    plot(t, P_selected(timeIdx, i), 'b-','LineWidth',1.5);
    plot([t(1), t(end)], [threshold(i), threshold(i)], 'r:','LineWidth',1.5);
    ylabel('P value');
    xlim([t(1), t(end)]);
    
    if i == 1
        legend({'Raw waveform', 'Selected waveform','Calculated P', 'Selected P', 'Threshold'}, ...
               'Location', 'north', 'Orientation', 'horizontal', ...
               'Position', [0.330041675855719,0.92639759908426,0.574999989196658,0.028749999292195]);
    end
    set(gca,'FontSize',12)
end



%% The benefits of waveform selection are compared by cross-correlation of these five stations

% the location of the five stations
lon = [-91.061897,-92.975899,-94.292999,-92.8535,-96.874496];
lat = [39.674999,39.045502,37.0354,36.3563,35.152699];

% the inter-station distance
R = [];   
RID = []; 
n=1;
for i = 1:nsta-1
    for j = i+1:nsta
        R(n) = distance(lat(i), lon(i), lat(j), lon(j), referenceEllipsoid('WGS 84'))./1000; % km
        RID(n,:) = [i, j];
        n=n+1;
    end
end
npairs=length(R);

% the number of sliding windows in the waveform data
win_len_samples = NW * stats.len_of_win;   
num_of_sld = stats.num_of_sld;


% calculate the cross-correlation of the raw waveform data
CC1 = 0;
CS1 = 0;
for i = 0:num_of_sld-1
    idx_start = i * win_len_samples + 1;
    idx_end   = (i+1) * win_len_samples;
    wave_win = waveform(idx_start:idx_end, :);
    [cc, cs] = cross_correlation(wave_win, RID);
    CC1 = CC1 + cc;
    CS1 = CS1 + cs;
end

% calculate the cross-correlation of the selected waveform data
CC2 = 0;
CS2 = 0;
for i = 0:num_of_sld-1
    idx_start = i * win_len_samples + 1;
    idx_end   = (i+1) * win_len_samples;
    wave_win = waveform_selected(idx_start:idx_end, :);
    wave_win(isnan(wave_win)) = 0;
    [cc, cs] = cross_correlation(wave_win, RID);
    CC2 = CC2 + cc;
    CS2 = CS2 + cs;
end

% plot the cross-correlation
Nt_win = win_len_samples; 
time_axis = (-Nt_win/2 : Nt_win/2 - 1) / FS; 

figure('Position', [100, 100, 1000, 400]);

subplot(1,2,1);
hold on;
for k = 1:npairs
    ydata = 2e12 * CC1(:,k) / CS1(k) + R(k);
    plot(time_axis, ydata, 'k');
end
xlim([time_axis(1) time_axis(end)])
xlabel('Time (s)');
ylabel('Inter-station distance (km)');
title('Cross-correlation of raw waveforms');

subplot(1,2,2);
hold on;
for k = 1:npairs
    ydata = 2e12 * CC2(:,k) / CS2(k) + R(k);
    plot(time_axis, ydata, 'k');
end
xlim([time_axis(1) time_axis(end)])
xlabel('Time (s)');
ylabel('Inter-station distance (km)');
title('Cross-correlation of selected waveforms');


% define a simple function to calculate the cross-correlation in frequency domain
function [cc, cs] = cross_correlation(wave, RID)
    Nt = size(wave, 1);
    Nf = floor(Nt/2) + 1;  
    F = fft(wave, [], 1); 
    F_pos = F(1:Nf, :);    
    
    npairs=length(RID(:,1));
    cc_full = zeros(Nt, npairs);
    for k = 1:size(RID,1)
        prod_pos = conj(F_pos(:, RID(k,1))) .* F_pos(:, RID(k,2));
        prod_full = [prod_pos; conj(prod_pos(end-1:-1:2, :))];
        c = ifft(prod_full, [], 1);
        cc_full(:, k) = ifftshift(c);
    end
    cs = sign(sum(abs(cc_full), 1));
    cc = cc_full;
end