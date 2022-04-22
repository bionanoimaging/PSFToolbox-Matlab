function xout = czt_z( xin , scalz, method )
%___________________________________________________________________________________
%
%  Aufruf:
%  xout = czt_two(xin,scalx,scaly,idir);
%
%  Beschreibung:
%  Berechnung der zweidimensionalen Chirp-z-Transformation mittels Faltung
%  mit Skalierungsfaktoren scalx bzw. scaly für ein zweidimensionales Feld
%  xin. Beide Koordinatenschnitte werden unabhängig behandelt. 
%  Das Feld kann rechteckig sein mit Dimensionen npx / npy.
%  Basis ist die Routine czt aus der Signal-Toolbox von Matlab.
%  Für scalx/scaly = 1 entsteht die FFT. Die Frequenzachsen werden um die Faktoren
%  scalx/scaly gespreizt.
%  Es muß scalx/scaly => 1 gelten.
%  Die Index-Shiftfunktion ist bereits enthalten, die Frequenz s = 0 liegt beim Index np/2+1.
%  Es gilt für die Frequenzrasterweite :    ds = 1 / ( npx * dx * scalx )  (y analog)
%
%  Version:
%  23.06.06   Herbert Gross  7.0  Überarbeitete Version 
%  28.11.18   Rainer Heintzmann and Rastimandresy Dina Miora
%  27.08.19   Rainer Heintzmann 
%
%  Input:   xin(npx,npy)  : 3D input field
%           scalz         : Scaling factor in Z
%           method : 'forward' or 'inverse'
%  Output:  xout : Output
%
%
%  Referenzen:
%  Rabiner, Schafer, Rader, The Cirp z-Transform Algorithm, IEEE Trans AU 17(1969) p. 86
%
%  Beispiel: Transformation eines Gaussprofils, verschiedene Punktzahlen und 
%  Spreizungsfaktoren
%
%  npx = 256  ; npy = 128 ;
%  scalx = 8 ;  scaly = 4 ; idir = 0 ;
%  a = 1; xmax = 2.;xmin = -xmax;
%  dx = ( xmax-xmin ) / (npx-1); dy = ( xmax-xmin ) / (npy-1);
%  for j=1:npx ;   x(j) = xmin + (j-1)*dx; end
%  for k=1:npy ;   y(k) = xmin + (k-1)*dy ; end  
%  [yp xp] = meshgrid(y,x); rp = sqrt(xp.^2+yp.^2);indin = find( rp < a );
%  xin = zeros(npx,npy,1); xin(indin) = 1 ;
%  xout = czt_2d( xin , scalx , scaly , idir );
%  figure(1);
%  pcolor(abs(xout));shading flat
%___________________________________________________________________________________
%

if nargin < 3
    method = 'forward';
end

if ~isa(xin,'double')
    xin = double(xin);
end

switch method 
    case'forward'
%        xin = xin;
    case 'inverse'
        xin = conj(xin);
    otherwise
        error('Unknown method. Choose ''forward'' or ''inverse'' as method.')
end

sz = size(xin);
npz = sz(3);
rampz = double(ramp([1,1,npz],3,'corner')); 
longrampz = double(ramp([1,1,npz*2-1],3)); % is also centerde 
longrampzone = double(ramp([1,1,npz*2-1],3,'corner'))+1; % is also centered 
% xout = zeros(npx,npy,1);
%
%  z-Richtung
%
%   fak = zeros(1,npz);
   n2 = npz/2+1;
   nn = rampz;
   kk = longrampz; % ( (-npz+1):(npz-1) ).';
   kk2 = (kk .^ 2) ./ 2;
%
   w = exp( -1i*2*pi/(npz*scalz));
   a = exp( -1i*pi/scalz );
   ww = w .^ kk2;
   aa = a .^ ( -nn );
   aa = aa .* ww(npz+nn);
   fv = fft( 1 ./ ww(longrampzone), 2*npz );   
   %fv = reshape(fv,[1,1,size(fv,1)]);
   %aa = reshape(aa,[1,1,size(aa,1)]);
%    for j=1:npx 
%        fak(j) =  exp(pi*1i*(j-n2)/scalz );
%    end
   fak =  exp(pi*1i*(rampz+1-n2)/scalz );
   fak = ww( npz+rampz) .* fak;
%   fak = reshape(fak,[1,1,size(fak,1)]);
%
   y = xin .* aa;
   g = ifft( fft(  y , 2*npz, 3) .* fv ,[],3);
   xout = ( g(:,:, npz:(2*npz-1) ) .* fak ) ;
%    for k=1:npy
%       y = xin(:,k) .* aa;
%       g = ifft( fft(  y , 2*npx ) .* fv );
%       xout(:,k) = ( g( npx:(2*npx-1) ) .* fak ) ;
%    end
%

if strcmp(method,'forward')
%    xout = xout; % Originally this is conj(xout) but I don't know yet why 
elseif strcmp(method, 'inverse')
%     xout = xout / ( npx*npy*scalx*scaly  );
    xout = conj(xout);
end

xout = dip_image(xout / sqrt(prod(sz)));   % modification to be compatible with the dipImage ft routine
