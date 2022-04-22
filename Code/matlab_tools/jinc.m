% res=jinc(mysize,myscales)  : Caculates a bessel(1,2*pi*radius)/radius   = jinc function, which describes the Airy pattern in scalar low NA approximation
%
% Example 1:
% pixelSize=203;  % half the Nyquist freq
% mysize=[256 256];
% lambda=488/pixelSize;
% na=0.3;
% AbbeLimit=lambda/na;   % coherent Abbe limit, central illumination, not incoherent
% ftradius=1/AbbeLimit*mysize;
% myscales=ftradius./mysize; % [100 100]/(488/0.3);
% res=jinc(mysize,myscales);  % Airy pattern (low NA) for these value (at n=1.0)
%
% Example 2: jinc such that the Fourier transformation has a defined radius myrad. E.g. needed for confocal pinholes
% mysize=[256 256];
% ftradius=10.0;  % pixels in Fourier space (or real space, if starting in Fourier space)
% myscales=ftradius./mysize; 
% ftpinhole=jinc(mysize,myscales);  % Fourier pattern for pinhole with radius ftradius
% ftpinhole=ftpinhole/ftpinhole(MidPosX(ftpinhole),MidPosY(ftpinhole))*pi*ftradius.^2;  % normalize
% pinhole=real(ft(ftpinhole))/sqrt(prod(mysize))  % normlized to one in the disk
function res=jinc(mysize,myscales)
if nargin < 1
    mysize=[256 256];
end
if nargin < 2
    pixelSize=203;  % half the Nyquist freq
    lambda=488/pixelSize;
    na=0.3;
    AbbeLimit=lambda/na;   % coherent Abbe limit, central illumination, not incoherent
    ftradius=1/AbbeLimit*mysize;
    myscales=ftradius./mysize; % [100 100]/(488/0.3);
end
myradius=pi*rrscale(mysize,myscales);
res=besselj(1,2*myradius) / (myradius);
res(MidPosX(res),MidPosY(res))=1;
