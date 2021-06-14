function IMGD=dilateN(IMG,MASK)
% IMG ... n-dim binary image (logical)
% MASK ... dilate n-dim mask (logical) or size of dilatation mask (pixels)
%
% Example:
% IMGD=dilateN(IMG,10) extends boundaries of 5 pixels 

if ~islogical(IMG)
   error('only for logical image') 
end    

if numel(MASK)==1
    MASK=single(ones(MASK*ones(1,length(size(IMG)))));
end

md=size(MASK);
if sum(mod(md,2)==0)>0 % eny dimension is odd
    warning('dilateN: mask has odd size')
end

IMG=single(IMG);
MASK=single(MASK);

try
    if gpuDeviceCount>0
%         disp('convn on GPU')
        IMGgpu=gpuArray(IMG);
        MASKgpu=gpuArray(MASK);
        IMGD=convn(IMGgpu,MASKgpu,'same');
        IMGD=gather(IMGD)>0;
    else
        IMGD=convn(IMG,MASK,'same')>0;
    end
catch
    IMGD=convn(IMG,MASK,'same')>0;
end

