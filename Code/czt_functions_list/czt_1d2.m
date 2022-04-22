% czt_1d(xin , scalx , scaled , d)
% 
% Chirp z transform along a single direction d of an ND array `xin` into the ND array 'xout'.
% Note that xin and xout can be the same array for inplace operations.
% This code is based on a 2D Matlab version of the CZT, written by H. Gross et al.
%     
% References: Rabiner, Schafer, Rader, The Cirp z-Transform Algorithm, IEEE Trans AU 17(1969) p. 86
% Modified by Dina R 22.08.2021
function xout = czt_1d2(xin, scaled, d,accuracy)

    sz=size(xin);
    dsize = sz(d);
    
    % accuracy non-zero should give better result but it seems not
    if nargin<4 || accuracy==0
        L = 2*dsize+2; % this value can be anything bigger than or equal to 2*dsize. Ideally it should be 2*ceil(log2(N+M-1)) to have fft most efficient but this can be costly in time, here we have N=M=dsize
    else
        L = 2^(ceil(log2(dsize)));
        fprintf(['Note: Higher accuracy with kernel of size equal to ',num2str(L),' is being used along dimension ',num2str(d),'\n']);
    end
    
    nn = (0:dsize-1); % vector with length equal to dsize
    AW = zeros(1,L);
    kk2 = AW; % length is equal to L 
    
    kk2(1:dsize) = nn.^2/2; 
    kk = L-dsize+1:L-1;
    kk2(kk+1) = ((nn(end:-1:2)) .^ 2) ./ 2; % vector same length as kk 
    
    w = exp(2*1i*pi/(dsize*scaled)) ; % scalar factor (tokony misy negative sign ve?)
    a = exp(1i*pi/scaled) ; % scalar factor (tokony misy negative sign ve?)
    ww = w .^ kk2; % vector same length as kk
    aa =  a.^ (-nn);
    AW(nn+1) = aa .* ww(nn+1) ; % is a 1d list of factors ; n = 0, 1 , ... , N-1
    V = fft(1./ww) ; % r = 0, 1, ..., L-1 // 1d list
    
    % instead of rotating the vectors to be along the desired dimension, set the desired dimension to be along the 1st dimension
    order = 1:length(sz);
    order( order == d ) = 1;
    order(1) = d;
    xinpermute = permute(xin,order); % change the orientation of the input image to match the dimension towards which we want the czt to run
    
    nsz = size(xinpermute); nsz(1) = L-dsize; % twice the size along direction d which is along the 1st dimension since we have rotate xin
    extraY=newim(nsz);
    xinpermute=cat(1,xinpermute,extraY);
    
    y = xinpermute*AW;
    ftboolean = zeros(1,length(nsz)); ftboolean(1) = 1;
    Y = dip_fouriertransform(y,'forward', ftboolean);  
    
    G = V*Y; 
    g = dip_fouriertransform(G,'inverse', ftboolean);

    roicenter = (sz - mod(sz,2))/2 ; % center of the region of interest
    xout = extract(g,sz,roicenter)*ww(nn+1);

    % xout is not centered and needs shift to the center
    szr = sz; szl = sz; 
    szl(1) = (dsize + mod(dsize,2))/2 ;  % size left
        cl = roicenter;
        cl(1) = (szl(1) + mod(szl(1),2))/2; % center of left 
    szr(1) = dsize - szl(1) ;  % size right
        cr = roicenter;
        cr(1) = cr(1) + (szr(1) + mod(szr(1),2))/2 ;  % center right

    xout = cat(1,extract(xout,szr,cr), extract(xout,szl,cl));
    
    xout = ipermute(xout,order); % permute data back to the original orientation
end


