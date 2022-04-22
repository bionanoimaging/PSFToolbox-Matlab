function imgout=ft2d(imgin)  % Two dimensional Fourier transform

transvec=zeros(1,ndims(imgin));
transvec(1:2)=1;

imgout=dip_fouriertransform(imgin,'forward',transvec);
