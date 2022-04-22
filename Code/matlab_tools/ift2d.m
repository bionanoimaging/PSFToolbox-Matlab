function imgout=ift2d(imgin)  % Two dimensional inverse Fourier transform

transvec=zeros(1,ndims(imgin));
transvec(1:2)=1;

imgout=dip_fouriertransform(imgin,'inverse',transvec);
