function [d,fs]=high_order_resample(d,fs,fs_new)


RFactor=fs/fs_new;
if RFactor>1
        
    decItterations=ceil(log10(RFactor));
    
    for i=1:decItterations
        if i==decItterations
            fs_out=fs_new;
        else
            fs_out=round((fs/RFactor^(1/decItterations)));
        end
        
%         sub_d_res=zeros(ceil(size(d,1)/(fs/fs_out)), size(d,2));
        ch1(:,1)=resample(d(:,1),fs_out,fs,100);
        sub_d_res=zeros(length(ch1), size(d,2));
        sub_d_res(:,1)=ch1; clear ch1;
        parfor ch=2:size(d,2) % døíve s parforem, nìkdy padá???
            sub_d_res(:,ch)=resample(d(:,ch),fs_out,fs,100);
        end
        d=sub_d_res; 
        fs=fs_out;
    end
else
   sub_d_res=[];
   parfor ch=1:size(d,2)
        sub_d_res(:,ch)=resample(d(:,ch),fs_new,fs,100);
   end
    d=sub_d_res;
    fs=fs_new;
end
