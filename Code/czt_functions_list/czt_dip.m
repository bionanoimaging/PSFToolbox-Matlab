function xout = czt_dip( xin , scalx , scaly , idir,method )
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
%
%  Input:   xin(npx,npy)  : Inputfeld, zweidimensional
%           scalx         : Lupenfaktor in x-Richtung
%           scaly         : Lupenfaktor in y-Richtung
%           idir          : Steuerparameter
%                           idir = 0 : Transformation in x und y
%                           idir = 1 : Transformation nur in x-Richtung
%                           idir = 2 : Transformation nur in y-Richtung
%
%  Output:  xout(npx,npy) : Outputfeld, zweidimensional
%
%  Abhängigkeiten : keine
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

if ~isa(xin,'double')
    xin = double(xin);
end

sz = size( xin );
if length(sz)>=3
    error('Choose a 2D input image')
end

switch method 
    case'forward'
        xin = xin;
    case 'inverse'
        xin = conj(xin);
    otherwise
        error("Unknown method. Choose 'forward' or 'inverse' as method.")
end

npx = sz(1);
npy = sz(2);
xout = zeros(npx,npy,1);
%
%  x-Richtung, 1. Index
%
if idir == 0 | idir == 1
%
   fak = zeros(1,npx);
   n2 = npx/2+1;
   nn = (0:(npx-1))';
   kk = ( (-npx+1):(npx-1) ).';
   kk2 = (kk .^ 2) ./ 2;
%
   w = exp( -i*2*pi/(npx*scalx));
   a = exp( -i*pi/scalx );
   ww = w .^ (kk2);   
   aa = a .^ ( -nn );
   aa = aa .* ww(npx+nn);
   fv = fft( 1 ./ ww(1:(2*npx-1)), 2*npx );   
   for j=1:npx ; fak(j) =  exp(pi*i*(j-n2)/scalx );end % size [1 256]
   fak = ww( npx:(2*npx-1) ) .* fak.' ; % size of >ww( npx:(2*npx-1) ) < is [256 1] and size of > fak.' < is [256 1] and size output is [256 1]
%
   for k=1:npy
      y = xin(:,k) .* aa; % same size // vector
      g = ifft( fft(  y , 2*npx ) .* fv );
      xout(:,k) = ( g( npx:(2*npx-1) ) .* fak ) ;
   end
%
else
    xout = xin ;
%
end
%
%  y-Richtung, 2. Index
%
if idir == 0 | idir == 2
%
   clear fak ;
   fak = zeros(1,npy);
   n2 = npy/2+1;
   nn = (0:(npy-1))';
   kk = ( (-npy+1):(npy-1) ).';
   kk2 = (kk .^ 2) ./ 2;
   w = exp( -i*2*pi/(npy*scaly));
   a = exp( -i*pi/scaly );
   ww = w .^ (kk2);   
   aa = a .^ ( -nn );
   aa = aa .* ww(npy+nn);
   fv = fft( 1 ./ ww(1:(2*npy-1)), 2*npy );   
   for j=1:npy ; fak(j) =  exp(pi*i*(j-n2)/scaly );end
   fak = ww( npy:(2*npy-1) ) .* fak.' ;
%
   for j=1:npx
      y = xout(j,:).' .* aa;
      g = ifft( fft(  y , 2*npy ) .* fv );
      xout(j,:) = ( g( npy:(2*npy-1) ) .* fak ).' ;
   end
%
end
%


if method == 'forward'
    xout = xout; % Originally this is conj(xout) but I don't know yet why 
elseif method == 'inverse'
%     xout = xout / ( npx*npy*scalx*scaly  );
    xout = conj(xout);
end
    
xout = dip_image(xout);
