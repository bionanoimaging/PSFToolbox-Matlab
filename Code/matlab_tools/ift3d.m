% res=ift3d(img)  : Performs a 3d inverse fft using the dip_fouriertransform function
function res=ift3d(img)  
if ndims(img) < 3
    error('ft3d needs at least a 3d input');
end
    resVec=zeros(1,ndims(img));
    resVec(1:3)=1;
    res=dip_fouriertransform(img,'inverse',resVec);
end
