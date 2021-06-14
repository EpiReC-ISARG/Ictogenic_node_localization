function [DTF,ttabs,settings]=avrDTF_MVAR_parametrisation_v01(d,tabs,fs,channel_err,bandwidth,ww,nn,mvar_max,ac_freq)
% DTF signal parameterisation. 
% The signal is sequentially filtered in defined bands and fragmented into
% time-segments. Average DTF outflow is computed for each time segment that is parallelized by the MATLAB pool (parpool).
% The final DTF matrix is assembled as products of the sub-bands. 
%   d ... iEEG signal (time x channels)
%   tabs ... time vector in datenum format (time x 1)
%   fs ... sampling frequency
%   channel_err ... vector of noise channels (1 x n)
%   bandwidth ... matrix of sub-bands to frequency fragmented DTF
%                 computation [f1low f1high; f2low f2high;...] f1<f2
%   ww ... window of time-segmentation in seconds.
%   nn ... overlap of time-segmentation in seconds, step=ww-nn
%   mvar_max... maximal MVAR order in bands. For narrow bands is used twice of bandwidth
%   ac_freq... frequency of AC-noise (e.g., 50 Hz - Europe, 60 Hz - USA)






% frequency band and MVAR order limits ====================================
bandwidth=round(bandwidth);
bandwidth(bandwidth<1)=1;
if max(bandwidth(:))>0.4*fs % cropping of bands over anti-aliasing filter below fs/2
    [row,~]=find(bandwidth<0.4*fs,1,'last');
    bandwidth=bandwidth(1:row,:);
end
disp(['fs=' num2str(fs) ' Hz, num of bands = ' num2str(size(bandwidth,1)) ', ' num2str(bandwidth(1,1)) '-' num2str(bandwidth(end,2)) ' Hz'])

% model order for each band
nAR=2*diff(bandwidth,1,2);
nAR=2*floor(nAR/2);% only odd (spectral simetry)
nAR(nAR>mvar_max)=mvar_max;



% signal pre-processing ===================================================
channel_err=sort(unique(channel_err));


d(:,channel_err)=[]; % removing of noise channels
md=size(d); % size of signal matrix

d=filt50Hz(d,fs,ac_freq,bandwidth(end,2)); % AC main hum filter (all mulltipliers by biquad notch) 50/60Hz



% time-segmentation settings ==============================================
winsize=round(ww*fs); % time-segmentation in samples
noverlap=round(nn*fs); % overlap in samples
index=1:winsize-noverlap:md(1)-winsize+1; % first indexes of time-segmentation  
winsizeT=ww; % in seconds
noverlapT=nn; % in seconds
index=index/fs; % in seconds



% CPUs parallelization ====================================================
% support a distributed computing (65 threads tested). The signal time-segments are divided into threads
try
    CPUpool=gcp; % get parpool (MATLAB Parallel Computing Toolbox)
    nCPU=CPUpool.NumWorkers;
catch
    nCPU=1; % without pool
end


if length(index)>nCPU % if numbers of segments is bigger than number of threads
    n_index=floor((length(index)/nCPU)); % number of segments per thread
    nCPU=floor(length(index)/n_index); 
else
    nCPU=length(index); % else use appropriate number of threads (for short signal or huge CPUs-pool)
end

% indexes of signal segments for individual threads
start_idx=zeros(1,nCPU);
start_seg=zeros(1,nCPU);
stop_idx=zeros(1,nCPU);
stop_seg=zeros(1,nCPU);
for i=1:nCPU
    start_idx(i) =  n_index*(i-1)+1;
    start_seg(i) = index(start_idx(i));
    
    if i==nCPU
        stop_idx(i) = length(index);
        stop_seg(i) =  md(1)/fs;
    else
        stop_idx(i) = n_index*i;
        stop_seg(i) =  index(stop_idx(i))+winsizeT-1/fs;
    end
end
% -------------------------------------------------------------------------


DTF    = single(zeros(md(2),max(bandwidth(:)),length(index))); % initialization of DTF matrix space 
ttabs=tabs(round(index*fs+winsize-1)); % datenum time of segments




% PROCESSING ==============================================================
seg_fs=fs; % variable of sampling frequency after signal downsampling
% DTF is sequentially computed for each band (from higher to lower frequency)
% Signal is sequentially downsampled for stable bandpass design and filtering
for i=size(bandwidth,1):-1:1 % DTF for each bands 
    
    % if higher-frequency of a band is higher than 1/3 of fs/2, signal is downsampled 
    if 3*bandwidth(i,2)<seg_fs/2
        r=seg_fs/(3*bandwidth(i,2));
        dec_fs=round(seg_fs/r);
    else
        dec_fs=seg_fs;
    end
    
    if seg_fs~=dec_fs
        d=high_order_resample(d,seg_fs,dec_fs); % resampling for factor also R>10
        seg_fs=dec_fs;
    end
    % --------------------------------------------------------------------- 
    
    
    
    
    tic;
    fprintf(1,['band  ' num2str(bandwidth(i,1)) '-' num2str(bandwidth(i,2)) ' Hz (' num2str(nCPU) ' seg.)... ']);
    % band-pass filter design (monotonic amplitude response - Butterworth) -----------------------------------------------------
    Wp=2*bandwidth(i,:)/seg_fs;
    Ws=2*[bandwidth(i,1)-0.1*seg_fs,  bandwidth(i,2)+0.1*seg_fs]/seg_fs; Ws(Ws<=0)=0.01; Ws(Ws>=1)=0.99;
    Rp=3;
    Rs=60;
    [n,Wn] = buttord(Wp,Ws,Rp,Rs);
    [bp,ap] = butter(n,Wn);
    
    % test of design - is monototic?
    [h,~]=freqz(bp,ap,1000,seg_fs);
    if sum(isnan(h))>0 || sum(abs(h)>(1+1e-6))
        error('filter design is wrong')
    end
    
    
    % Preparing of data parcelation for each thread d->SEG{} +2s around to removing impulse responses of filter (at start and end)
    % Using of SEG{} array reduces transferred data to CPUs-pool through LAN. Only appropriate iEEG epoch is transfered.  
    for s=1:nCPU
        if s==1
            SEG{s}=d(ceil(start_seg(s)*seg_fs):round(stop_seg(s)*seg_fs+2*seg_fs),:);
        elseif s==nCPU
            SEG{s}=d(ceil(start_seg(s)*seg_fs-2*seg_fs):round(stop_seg(s)*seg_fs),:);
        else
            SEG{s}=d(ceil(start_seg(s)*seg_fs-2*seg_fs):round(stop_seg(s)*seg_fs+2*seg_fs),:);
        end
    end
   
    
    % =====================================================================
    % Computation of average DTF outflow from channel j->all for the band, for each time-segment (ww) with overlap (nn)
    % Computation is parallelized in parpool for time-segments (parfor) 
    parfor s=1:nCPU % for CPU thread 
        % band pass filtering ---------------------------------------------
        if s==1 % for section at the beginning of iEEG
            SIG=SEG{s};
            SIG=filtfilt(bp,ap,SIG); % zero-phase filtering
            SIG=SIG(1:end-2*seg_fs,:); % removing of impulse responses on end
            
        elseif s==nCPU % for section at the middle of iEEG
            SIG=SEG{s};
            SIG=filtfilt(bp,ap,SIG); % zero-phase filtering
            SIG=SIG(2*seg_fs+1:end,:); % removing of impulse responses
            
        else
            SIG=SEG{s}; % for section on end of iEEG
            SIG=filtfilt(bp,ap,SIG); % zero-phase filtering
            SIG=SIG(2*seg_fs+1:end-2*seg_fs,:); % removing of impulse responses in the begining
        end
        
        
        Nsub=stop_idx(s)-start_idx(s)+1;
        subindex=round([1 (1:Nsub-1)*(winsizeT-noverlapT)*seg_fs]); % time-indexes corresponding to signal
        
        % MAIN DTF-parametrization ========================================
        s_DTF{s}=direct_function(SIG,seg_fs,[bandwidth(i,1), bandwidth(i,2)],floor(winsizeT*seg_fs),floor(noverlapT*seg_fs),nAR(i),subindex);
        % =================================================================
    end
    
    
    
    % Serialization of threads
    for s=1:nCPU
        DTF(:,bandwidth(i,1):bandwidth(i,2),start_idx(s):stop_idx(s))=s_DTF{s};
    end
    
    Ttoc=toc;
    
    if Ttoc<60
        fprintf(1,['pass after ' num2str(Ttoc,'%.1f') ' seconds\n']);
    elseif Ttoc<3600
        fprintf(1,['pass after ' num2str(Ttoc/60,'%.1f') ' minutes\n']);
    else
        fprintf(1,['pass after ' num2str(Ttoc/3600,'%.1f') ' hourss\n']);
    end
end   


% Extension of DTF matrix (ch x band x time) to add the initially excluded noise channels
% the noise channels are replaced by zeros
for ch=channel_err
    DTF=[DTF(1:ch-1,:,:);zeros(1,size(DTF,2),size(DTF,3));DTF(ch:end,:,:)];
end

% Stack vector corresponding with analyzed frequencies (1) and skipped frequencies (0)
F=false(max(bandwidth(:)),1);
for i=1:size(bandwidth,1)
    F(bandwidth(i,1):bandwidth(i,2))=true; % indicator of computed frequency line
end

% Structure of used settings
settings.winsize=ww;
settings.noverlap=nn;
settings.fs=fs;
settings.nAR=nAR;
settings.channel_err=channel_err;
settings.bandwidth=bandwidth;
settings.ac_freq=ac_freq;
settings.F=F;

end











function DTF=direct_function(SIG,fs,bandwidth,winsize,noverlap,AIC_r,index)
% MVAR model estimation and DTF computation
warning off all



DTF   =single(zeros(size(SIG,2),diff(bandwidth)+1,length(index))); % memory reservation
f=bandwidth(1):1:bandwidth(end); % frequency vector of the band 

NAN_count=0; % counter of MVAR estimation failure
for k=1:length(index)
    i=index(k);
    
    segment=SIG(i:i+winsize-1,:); % iEEG time-segment

    
    
    [AR,~,PE]=mvar_simplified(segment,AIC_r); % MVAR model estimation ===================
    if isnan(AR(1,1))
        NAN_count=NAN_count+1;
        continue 
    end
     
    % MVAR -> DTF =========================================================
    [~, ~, ~, ~, full_DTF] = mvfreqz_faster(eye(size(segment,2)),[eye(size(segment,2)),-AR],PE(:,AIC_r*size(segment,2)+(1:size(segment,2))),f,fs);
    
    
    % Average outflow j->all ---------------------------------------------
    DTF(:,:,k)  =permute(single(sum(full_DTF.*(1-eye(size(full_DTF,1))),1)/(size(full_DTF,1)-1)'),[2 3 1]); % (ch x freq. x time)


end

if NAN_count>0
    disp(['MVAR model failed - ' num2str(NAN_count) 'X (' num2str(((NAN_count-1)*(winsize-noverlap)+winsize)/fs) 'sec.)']);
    disp('Check iEEG signal and remove the correlated channels (e.g. Corr>0.99)')
end
end
    