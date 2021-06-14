function d=filt50Hz(d,fs,hum_f,f_max)

if ~(exist('hum_fs')==1);    hum_f=50; end
if isempty(hum_f); hum_f=50; end


if ~(exist('f_max')==1)
    f_max=fs/2;
end

% if ~(exist('M')==1)
%     M=matlabpool('size'); if M==0; M=1; end
% end



if min(size(d))==1
   d=d(:); 
end

R = 1; r = 0.985;


f0 = hum_f:hum_f:f_max; % Hz


for i=1:length(f0)
    b = [1 -2*R*cos(2*pi*f0(i)/fs) R*R];
    a = [1 -2*r*cos(2*pi*f0(i)/fs) r*r];
%     parfor (ch=1:size(d,2),M)
%         d(:,ch)=filtfilt(b,a,d(:,ch));
%     end
   d=filtfilt(b,a,d);
end


