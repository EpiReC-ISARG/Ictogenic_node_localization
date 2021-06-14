function IMGE=erodeN(IMG,MASK)
% IMG ... n-dim binary image (logical)
% MASK ... erode n-dim mask (logical) or size of dilatation mask (pixels)
%
% Example:
% IMGE=dilateN(IMG,10) crops boundaries of 5 pixels 

if ~islogical(IMG)
   error('only for logical image') 
end    

if numel(MASK)==1
    MASK=true(MASK*ones(1,length(size(IMG))));
end

md=size(MASK);
if sum(mod(md,2)==0)>0 % eny dimension is odd
    warning('erodeN: mask has odd size')
end


IMG=single(IMG);
MASK=single(MASK);

try
    if gpuDeviceCount>0
%         disp('convn on GPU')
        NIMG=gpuArray(abs(IMG-1));
        MASKgpu=gpuArray(MASK);
        IMGE=convn(NIMG,MASKgpu,'same');
        IMGE=~logical(gather(IMGE));
    else
        warning('GPU error')
        IMGE=~logical(convn(~IMG,MASK,'same'));
    end
catch
    IMGE=~logical(convn(~IMG,MASK,'same'));
end
