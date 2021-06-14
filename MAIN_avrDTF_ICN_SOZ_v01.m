clc; close all; clear all;

%% iEEG data loading
data=load('P36_2011-06-19_17-46_seizure.mat');
noisy_channels=[17 42 55 58 59 64 65]; % artefacts channels after visual inspection
ac_freq=50; % main hum frequency (Hz)

addpath('mfiles')
data.d=double(data.d).*repmat(double(data.mults),[size(data.d,1) 1]); % iEEG in double (time x ch)
% NOTE:
% data.fs; % sampling frequency
% data.tabs; % time vector in datenum format (in column: time x 1)



%% avrDTF parametrization
% simplified algorithm for ictogenic nodes (ICNs) localization

bandwidth=[2 12;12 25;25 48; 52 75; 75 98; 102 148;152 198;202 248; 252 298]; % freq. bands for DTF distributed computing
% NOTE: The algorithm contains AC-hum filters, but we recommend using "bandwidth"
% skipping multiple of main hum frequency, e.g. for 50 Hz: [...48; 52... ; ...98; 102...]    


ww=1; % (in sec.) 1 second segmentation 
nn=0.9; % (in sec.) overlap 0.9 seconds (0.1 s steps)
mvar_max=20; % MVAR orders are twice of bandwidth or maximally 20. 



% Warning: DTF parametrisation can be high CPU and memory demanding. 
% 1) We recomend >1GB RAM per CPU-core/thread. 
% 2) If number of channels >>100, you consider channels reduction or increase segmentation window (ww=1.5-2 s), whorse time-resolution 
% NOTE: The parameterisation supports parallel distributed computing by parpool (tested in MATLAB 2020a, 65 threads)
[DTF,ttabs,settings]=avrDTF_MVAR_parametrisation_v01(data.d,data.tabs,data.fs,noisy_channels,bandwidth,ww,nn,mvar_max,ac_freq); % MAIN DTF FUNCTION


% saving of DTF parameterisation
save(['avrDTF_' datestr(data.tabs(1),'yyyymmdd-HHMMSS') '.mat'],'DTF','ttabs','settings')


%% seizure probability and ICNs detection for prospective assesment
% load(['avrDTF_' datestr(data.tabs(1),'yyyymmdd-HHMMSS') '.mat'])

soz=[9 31 32 26 41 53 56 57]; % channels in seizure onset zones marked by clinicians (not required)
sot=data.sot; % seizure onset time (by clinicians) in datenum format. The accurate time is not necessary, but the value determines end of pre-ictal period for training. 
% In practice, the time can be set 10-30 s before seizure semiology to define non-ictal activity. 
% The dramatic increase in seizure probability indicates seizure onset time without clinical knowledge, see figure(2), subplot(2,2,3). 
ma_window=3; % moving average window (in seconds) to compute a prior probability of seizure-related ICNs
decreasing_band=[25 175]; % (Hz) frequency range of specific DTF decreasing for ICNs detection
th=0.25; % detecting threshold
nu=0.05; % 5th percentile of outlier (non-ictal variability is between 5-95%, ictal is outside)


% signal vizualization ====================================================
figure(1); clf; % iEEG  
% signal filtering --------------------------------------------------------
data.d=data.d-repmat(trimmean(data.d,75,2),[1 size(data.d,2)]); % removing of average reference componets
[bhp,ahp]=butter(5,2*2/data.fs,'high'); % IIR hig-pass filter >2Hz
data.d=filtfilt(bhp,ahp,data.d); % zero phase filtering
data.d=filt50Hz(data.d,data.fs,ac_freq,data.fs/2); % AC noise filtering
data.d(:,noisy_channels)=0; % hide/zeros noise channels
% multichannel plotting ---------------------------------------------------
DEL=0.5*median(max(abs(data.d(:)))); % channel distancing
YA=linspace(DEL*(size(data.d,2)-1),0,size(data.d,2)); % channel offset
p1=plot((data.tabs-sot)/datenum(0,0,0,0,0,1),data.d+repmat(YA,size(data.d,1),1),'k'); % iEEG
axis tight; title('iEEG signal with clinical assesment'); xlabel('time to seizure (s)')
set(gca,'YTick',YA(end:-1:1),'YTickLabel',data.header.label(end:-1:1))
hold on
p2=plot([0 0],YA([1 end]),':r');
p3=plot(0*ones(length(soz),1),YA(soz),'r.','MarkerSize',15);
legend([p1(1) p2 p3],{'iEEG (>2 Hz)','seizure onset','SOZ'})





% Identification of ICNs ==================================================

% determination of pre-ictal variability (non-epileptic connectivity) 
x_baseline=ttabs<sot; % pre-ictal segments
ql_DTF=quantile(DTF(:,:,x_baseline),nu,3); % lower outliers threshold q0.05

% epiletiform activity-related outliers in DTF matrix (true/false)
c_decrease= DTF < repmat(ql_DTF,[1 1 size(DTF,3)]); % decrease in  connectivity

% a prioro probability in moving widow (3 s)
maorder=round((ma_window-settings.winsize)/(settings.winsize-settings.noverlap))+1; % order of a prior probability MA-filter
bb=ones(1,1,maorder); bb=bb./numel(bb); % (1 x 1 x time) - coefficients of MA-filter
apriory_decrease=convn(c_decrease,bb,'same'); % a prior probability by convolution (FIR: MA-filter)

% average probability of DTF decrease in specific frequency range [25 175] to determone a seizure probability 
ICN_probability=sum(apriory_decrease(:,decreasing_band(1):decreasing_band(2),:),2)/sum(settings.F(decreasing_band(1):decreasing_band(2)));
ICN_probability=permute(ICN_probability,[1 3 2]);


time_rel=(ttabs-sot)/datenum(0,0,0,0,0,1); % relative time-axis (sot = 0s)


% Visualization of seizure probability ====================================
figure(2); clf; 
subplot(221) 
imagesc(time_rel,1:size(data.d,2),ICN_probability); % a probability of DTF decreasing (ch x time)
cb=colorbar; cb.Label.String='probability of DTF decrease';
title(['a probability of seizure in ' num2str(decreasing_band(1)) '-' num2str(decreasing_band(2)) 'Hz'])
xlabel('time to seizure (s)')
set(gca,'YTick',1:size(data.d,2),'YTickLabel',data.header.label(1:end),'fontsize',8)
colormap('jet')

subplot(222) % (true/false) a probability > ICNs detecting threshold (ch x time)
imagesc(time_rel,1:size(data.d,2),ICN_probability>th); 
title(['a probability >' num2str(100*th) '%'])
xlabel('time to seizure (s)')
set(gca,'YTick',1:size(data.d,2),'YTickLabel',data.header.label(1:end),'fontsize',8)


% automatic seizure detection =======================================================
pp_max=max(ICN_probability,[],1); % maximal ICNs probability over channels to identify seizure onset time
subplot(223)
plot(time_rel,pp_max); axis tight; ylim([0 1]); title('maximal seizure probability over channels')
hold on
plot(time_rel([1 end]),[th th],'r')
xlabel('time to seizure (s)')

% detection of seizure onset ----------------------------------------------
% The seizure is characterised by dramatic increase in probability. Suspect epochs <3s can be removed
% by morphological operations (erosion and dilatation)  
pp_max_eroded=erodeN(pp_max>0.25,round(3/(settings.winsize-settings.noverlap))/2); % removing of short epochs
pp_max_dilated=dilateN(pp_max_eroded,round(3/(settings.winsize-settings.noverlap))/2); % restoring of long epochs


% seizure onset time by connectivity change *******************************
sot_by_ICNs=find(pp_max_dilated,1,'first'); % find of 1st threshold crossing
% if maximal probability is pre-ictally varivale due to interictal
% discharges etc., you can select time to determine ICNs manually by the
% mouse (uncomment below ginput):

% [~,manual_soz_time]=ginput(1); sot_by_ICNs=find(time_rel>=manual_soz_time,1);

% *************************************************************************

disp(['seizure onset was detected at:' datestr(ttabs(sot_by_ICNs))])


plot(time_rel(pp_max_dilated),pp_max_dilated(pp_max_dilated>0),'k.'); ylim([0 1.05])
plot(time_rel(sot_by_ICNs)*[1 1],[0 1],':k')
legend('maximal probability','threshold','seizure','detected seizure onset','location','best')


% Identification of ICNs ================================================== 
subplot(222); hold on
rectangle('Position',[time_rel(sot_by_ICNs) 0.5 6.5 size(data.d,2)+0.5],'EdgeColor','r') % area of ICNs detection
text(time_rel(sot_by_ICNs),5,'area of ICN detection','color','r','HorizontalAlignment','right')

% DTF decresing with probability >25% at time up to 6.5s after the first detectable decrease in connectivity 
ICNs=find(sum(ICN_probability(:,time_rel>=time_rel(sot_by_ICNs) & time_rel<=(6.5+time_rel(sot_by_ICNs)))>th,2)>0); % ICNs

disp(['ICNs channel:' num2str(ICNs')])
disp(['ICNs channel:' cell2mat(arrayfun(@(x) [cell2mat(data.header.label(x)) ', '], ICNs,'UniformOutput',0)')])


% showing of clinically marked SOZ and detected ICNs channel  
subplot(222)
hold on
plot(time_rel(sot_by_ICNs)*ones(length(ICNs),1),ICNs,'c.','MarkerSize',15);
plot(0*ones(length(soz),1),soz,'r.','MarkerSize',15);
xlabel('time to seizure (s)')
legend('ICNs','SOZ','location','northwest')


% showing of relative connectogram in channels detected as ICNs
DTF_relative=DTF./repmat(mean(DTF(:,:,x_baseline),3),[1 1 size(DTF,3)]); % normalization by pre-ictal setting

ICN_connectogram=mean(DTF_relative(ICNs,:,:),1); % mean over ICN channels
ICN_connectogram=permute(ICN_connectogram,[2 3 1]); 
ICN_connectogram(isnan(ICN_connectogram))=1; % main AC hum (n x 50Hz) NaN->1

subplot(224)
imagesc(time_rel,1:size(ICN_connectogram,1),ICN_connectogram); caxis([0 2]); % the connectogram ICNs-outflow
axis xy; ylabel('f(Hz)'); title('average normalized connectogram of detected ICNs')
cb=colorbar; cb=colorbar; cb.Label.String='relative change in DTF';
cb.TickLabels{1}=[cb.TickLabels{1} ' - decrease'];
cb.TickLabels{end}=['>' cb.TickLabels{end} ' - increase'];


% add markers of ICNs to iEEG signal
figure(1);
hold on; 
p4=plot(time_rel(sot_by_ICNs)*[1 1],YA([1 end]),':c');
p5=plot(time_rel(sot_by_ICNs)*ones(length(ICNs),1),YA(ICNs),'c.','MarkerSize',15);
legend([p1(1) p2 p3 p4 p5],{'iEEG (>2 Hz)','seizure onset','SOZ','detected onset','ICNs'})