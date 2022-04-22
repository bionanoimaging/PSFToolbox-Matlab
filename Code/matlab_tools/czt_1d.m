% xout =  czt_1d(xin , scaled , d)
% 
% Chirp z transform along a single direction d of an ND array `xin` into the ND array 'xout'.
% scaled is a zoom-in or zoom-out factor.
% Note that xin and xout can be the same array for inplace operations.
% This code is based on a 2D Matlab version of the CZT, written by H. Gross et al.
%     
% References: Rabiner, Schafer, Rader, The Chirp z-Transform Algorithm, IEEE Trans AU 17(1969) p. 86
%
% Example:
%         xin = readim;
%         xout = czt_1d(xin, 2, 1);
%         xin0 = extract(xin,2*size(xin));
%         fxout = extract(ft1d(xin0),size(xin));
%         cat(3,xout,fxout) % czt zoomed twice , ft from zero-padded
% 
% Last edit: 29.08.21 R. Dina and R. Heintzmann

function xout = czt_1d(xin, scaled, d)

    sz=size(xin);
    % rtype = real(eltype(xin));% what is this doing?
    dsize = sz(d);
    nn = (0:dsize-1); % vector with length equal to dsize
    kk = ((-dsize+1):(dsize-1)); % length is equal to 2*dsize-1 // length is odd and vector is centered at dsize (zero-position)
    kk2 = (kk .^ 2) ./ 2; % vector same length as kk
    
    w = exp(-2*1i*pi/(dsize*scaled)) ; % scalar factor
    half_pix_shift = 2*floor(dsize/2)/dsize;
    a = exp(-1i*pi/scaled*half_pix_shift) ; % scalar factor
    ww = w .^ kk2; % vector same length as kk
    aa =  a.^ (-nn);
    aa = aa .* ww(dsize + nn) ; % is a 1d list of factors
    tofft = 1./ww(1:2*dsize-1);

    to_fft = extract(dip_image(tofft),2.*dsize,dsize);
    fv = ft(to_fft);% is always 1d

    y = xin*reorient(aa,d);
    nsz = sz; nsz(d) = 2*sz(d); % twice the size along direction d which is along the 1st dimension since we have rotate xin
    
    to_fft = extract(y,nsz,floor(nsz/2)); 
    ftboolean=zeros(1,length(nsz)); ftboolean(d)=1;
    g = dip_fouriertransform(dip_fouriertransform(to_fft,'forward',ftboolean).*reorient(fv,d), 'inverse', ftboolean);
    
    oldctr = floor(sz(d)/2)- 1;  % ADD -1

    if mod(dsize,2) == 1 % This is to deal with a strange phase shift appearing for odd-sized arrays
        extra_phase = (2*dsize - 2)/(2*dsize); % # 5: 12 / 15, 7: 12/14, 9: 16/18, 11: 20/22
    else 
        extra_phase = 1;
    end
    
    fak = ww(dsize:(2*dsize - 1)).* exp(1i*pi*xx(dsize,1)/scaled * extra_phase ); % is a 1d list of factors
    
    xout = extract(g,sz,oldctr).*reorient(fak,d);
   
    
    if mod(dsize,2) == 0 && scaled>1
         % set a vector positions
        invec=zeros(1,length(sz));
        for k=1:length(sz)
            if k==d
                invec(k)=1;
            else
                invec(k)=':';
            end
        end
        
        midp = floor(dsize/2); 
        for o = (1 + mod(midp,2):2:dsize)-1
            invec(d)=o;
            xout(invec) = xout(invec) - slice(xin,d, 1) .* (1i).^mod(o-midp,4);
        end
    end
    
end


