% [de,noiseFactor,offsetFactor]=DampEdge (im,part,dim,shape,method, sigma) : dims edges of the input image - Hann window with with central flat are
%
% INPUT:
% im - image which edges will be dimmed down
% part - percentage of the (shorter) side - this will be size of the are that will be dimmed down by windowing
% dim - number of dimensions that will be dimmed (can be 1, 2 or 3). When using a 4D image, each 3D volume can be dimmed, but not along the 4th dimension.
% shape - shape of window-function: shape=0 round shape,shape=1-square shape (default = 1) 
% method - 1: dim down to zero, 2: fade to a Gaussian blurred average of the opposing lines, 3: fade to a Gaussian blurred average but using a cos rather than a cos^2 function
% sigma : Optional size of the Gaussian blurring kernel for method 2, (default: 10)
%
% OUTPUT:
% de : result with edges damped
% noiseFactor : yields the fraction of remaining noise variance , if I had no offset
% offsetFactor : yields the fraction of offset-only remaining noise variance
%
% Example:
% damped = DampEdge(readim,0.05,2,1,2)

function [de,noiseFactor,offsetFactor]=DampEdge (im,part,dim,shape,method,sigma)
if nargin<2 || isempty(part)
    part=0.15; % dim at 15%
end
if nargin<3 || isempty(dim)
    dim=max(ndims(im),2);  % by default use 2D operation on 3D image
end
    
part=part*2; %this is for half of the side...
origdims=ndims(im);
mypref=dipgetpref('DerivativeFlavour');

if nargin<4 || isempty(shape)
    shape = 1; %     'square shape';
end
if nargin<5 || isempty(method)
    method = 3; % damp to average
end
if nargin<6 || isempty(sigma)
    sigma = 10;  % for Gaussian blurring
end

if (part == 0)
    part=1e-10;
end

if ~isreal(im)
    if nargout > 1
        [der,noiseFactor]= DampEdge (real(im),part,dim,shape,method,sigma);
        [dei,noiseFactor]= DampEdge (imag(im),part,dim,shape,method,sigma);
    else
        der = DampEdge (real(im),part,dim,shape,method,sigma);
        dei = DampEdge (imag(im),part,dim,shape,method,sigma);
    end
    de=der+i*dei;
    return;
end

sizevector = size(im);
sizevector3D = size(im);
sizevector2D(1) = sizevector(1);
if length(sizevector2D)> 1
    sizevector2D(2) = sizevector(2);
end

centervector=floor(sizevector/2);

sigmaX = sigma;
sigmaY = sigma;
sigmaZ = sigma;
if dim > ndims(im)
    dim = ndims(im);
end

if dim < 1 || size(im,1) < 2
    sigmaX=0;
end
if dim < 2 || size(im,2) < 2
    sigmaY=0;
end
if dim < 3 || size(im,3) < 2
    sigmaZ=0;
end

if method == 1  % dim down to zero
    if length(sizevector3D)>2 && dim==2
        im=squeeze(im);
        sizevector=size(im);
    end

    if shape == 0; %round shape
        out = min(sizevector2D);
        in = out * (1 - part);

        if dim > 2
            de1 = rr([sizevector],'freq');
        elseif dim == 2
            de1 = rr([sizevector(1:2)],'freq');
        else
            de1 = rr([sizevector(1)],'freq');
        end
        if numel(sizevector) > dim 
           de1=repmat(de1,[ones(1,dim) sizevector(dim+1:end)]);
        end
        mask = de1 < 0.5;
        de1 = (de1 - (in/out * 0.5)) * pi/0.5 * out/(out-in);
        de1(de1<0) = 0;
        de = (1 - 0.5 * (1 - cos(de1))) .* mask;

    elseif shape == 1 %square shape
        out = min(sizevector2D);
        in = out * (1 - part);
        l = out - in;

        out = sizevector(1);
        in = out - l;
        de1 = (abs(xx([sizevector],'freq')) - (in/out * 0.5)) * pi/0.5 * out/(out-in);
        de1(de1 <0) = 0;

        out = sizevector(2);
        in = out - l;
        de2 = (abs(yy([sizevector],'freq')) - (in/out * 0.5)) * pi/0.5 * out/(out-in);
        de2(de2 <0) = 0;
        
        if (dim > 2)
            out = sizevector(3);
            in = out - l;
            de3 = (abs(zz([sizevector],'freq')) - (in/out * 0.5)) * pi/0.5 * out/(out-in);
            de3(de3 <0) = 0;
        else
            de3=0;
        end

        mask = de1 - de2;
        mask = mask > 0;
        de1 = de1 * mask;
        de2 = de2 * ~mask;

        mask = de1 - de3;
        mask = mask > 0;
        de1 = de1 * mask;
        de3 = de3 * ~mask;

        mask = de2 - de3;
        mask = mask > 0;
        de2 = de2 * mask;
        de3 = de3 * ~mask;

        de = pi - (de1 + de2 + de3);
        de = 1  - (0.5 * (1 + cos(de)));
    else
        error('Shape can only me round (0) or rectanglular (1)');
    end
    de = im .* de;
     %de = im .* de;    
elseif method >= 2 
    dipsetpref('DerivativeFlavour','fourier');  % ensure the gaussian filtering is a cyclic convolution
    if shape == 0
        error('Gaussian blur of the edges is only allowed for square shape (shape=1)');
    end
    if ndims(im) < 4 && size(im,3) <= 1
       im = squeeze(im);
    end
    im2 = im*0;
    sumw=im2+0;
    maxw=im2+0;
    if ndims(im) == 3 && size(im,3) > 1
        % error('3d dimming not implemented yet');
        if (dim > 1)
            border=ceil(size(im,1)/2*(part/2));
            if (border > 1)
                if (method>=3)
                    borderdivisor=(border-0.5);  % This is used for the right scaling of the weights near the edges
                else
                    borderdivisor=(border-1);  % For compatibility reasons this is kept (method 2). This means that the left and right pixels are both identical
                end
            else
                borderdivisor=1.0;
            end
            w=xx([border,size(im,2),size(im,3)],'corner')/borderdivisor;  % make a linear weigth from left to right
            if (method == 3)
                w=(1-cos(w*pi/2));  % cos ^2 weighting
            else
                w=(1-cos(w*pi))/2;  % cos ^2 weighting
            end
            wi=w(end:-1:0,:,:);
            % Avg=repmat(permute(expanddim(gaussf(squeeze(im(0,:,:)/2+im(end,:,:)/2),[sigmaY,sigmaZ]),3),[3 1 2]),[border 1 1]);
            Avg=permute(expanddim(gaussf(squeeze(im(0,:,:)/2+im(end,:,:)/2),[sigmaY,sigmaZ]),3),[3 1 2]);
            tmp=Avg .* wi + im2(0:border-1,:,:); 
            im2(0:border-1,:,:)=tmp;
            tmp=sumw(0:border-1,:,:)+wi;
            sumw(0:border-1,:,:)=tmp;
            tmp=max(maxw(0:border-1,:,:),wi);
            maxw(0:border-1,:,:)=tmp;
            tmp=Avg .* w+im2(end-(border-1):end,:,:);
            im2(end-(border-1):end,:,:)=tmp;
            tmp=sumw(end-(border-1):end,:,:)+w;
            sumw(end-(border-1):end,:,:)=tmp;
            tmp=max(maxw(end-(border-1):end,:,:),w);
            maxw(end-(border-1):end,:,:)=tmp;
        end
        if (dim >= 2)
            border=ceil(size(im,2)/2*(part/2));
            if (border > 1)
                if (method>=3)
                    borderdivisor=(border-0.5);  % This is used for the right scaling of the weights near the edges
                else
                    borderdivisor=(border-1);  % For compatibility reasons this is kept (method 2). This means that the left and right pixels are both identical
                end
            else
                borderdivisor=1.0;
            end
            w=yy([size(im,1),border,size(im,3)],'corner')/borderdivisor;  % make a linear weigth from left to right
            if (method == 3)
                w=(1-cos(w*pi/2));  % cos ^2 weighting
            else
                w=(1-cos(w*pi))/2;  % cos ^2 weighting
            end
            wi=w(:,end:-1:0,:);
            % Avg=repmat(permute(expanddim(gaussf(squeeze(im(:,0,:)/2+im(:,end,:)/2),[sigmaX,sigmaZ]),3),[1 3 2]),[1 border 1]);
            Avg=permute(expanddim(gaussf(squeeze(im(:,0,:)/2+im(:,end,:)/2),[sigmaX,sigmaZ]),3),[1 3 2]);
            tmp=im2(:,0:border-1,:)+Avg.*wi; 
            im2(:,0:border-1,:)=tmp;
            tmp=sumw(:,0:border-1,:)+wi;
            sumw(:,0:border-1,:)=tmp;
            tmp=max(maxw(:,0:border-1,:),wi);
            maxw(:,0:border-1,:)=tmp;
            tmp=im2(:,end-(border-1):end,:)+Avg.*w;
            im2(:,end-(border-1):end,:)=tmp;
            tmp=sumw(:,end-(border-1):end,:)+w;
            sumw(:,end-(border-1):end,:)=tmp;
            tmp=max(maxw(:,end-(border-1):end,:),w);
            maxw(:,end-(border-1):end,:)=tmp;
        end
        if (dim >= 3) && (size(im,3) > 1)
            border=ceil(size(im,3)/2*(part/2));
            if (border > 1)
                if (method>=3)
                    borderdivisor=(border-0.5);  % This is used for the right scaling of the weights near the edges
                else
                    borderdivisor=(border-1);  % For compatibility reasons this is kept (method 2). This means that the left and right pixels are both identical
                end
            else
                borderdivisor=1.0;
            end
            w=zz([size(im,1),size(im,2),border],'corner')/borderdivisor;  % make a linear weigth from left to right
            if (method == 3)
                w=(1-cos(w*pi/2));  % cos ^2 weighting
            else
                w=(1-cos(w*pi))/2;  % cos ^2 weighting
            end
            wi=w(:,:,end:-1:0);
            % Avg=repmat(permute(expanddim(gaussf(squeeze(im(:,:,0)/2+im(:,:,end)/2),[sigmaX,sigmaY]),3),[1 2 3]),[1 1 border]);
            Avg=permute(expanddim(gaussf(squeeze(im(:,:,0)/2+im(:,:,end)/2),[sigmaX,sigmaY]),3),[1 2 3]);
            tmp=im2(:,:,0:border-1)+Avg.*wi; 
            im2(:,:,0:border-1)=tmp;
            tmp=sumw(:,:,0:border-1)+wi;
            sumw(:,:,0:border-1)=tmp;
            tmp=max(maxw(:,:,0:border-1),wi);
            maxw(:,:,0:border-1)=tmp;
            tmp=im2(:,:,end-(border-1):end)+Avg.*w;
            im2(:,:,end-(border-1):end)=tmp;
            tmp=sumw(:,:,end-(border-1):end)+w;
            sumw(:,:,end-(border-1):end)=tmp;
            tmp=max(maxw(:,:,end-(border-1):end),w);
            maxw(:,:,end-(border-1):end)=tmp;
        end
        tmp=sumw>0;
        im2(tmp)=im2(tmp)./sumw(tmp);
        de=im2*maxw+im*(1-maxw);
    elseif ndims(im)<=3 
        if (dim == 1 && ndims(im) < 2)
            border=ceil(size(im,1)/2*(part/2));
            if (border > 1)
                if (method>=3)
                    borderdivisor=(border-0.5);  % This is used for the right scaling of the weights near the edges
                else
                    borderdivisor=(border-1);  % For compatibility reasons this is kept (method 2). This means that the left and right pixels are both identical
                end
            else
                borderdivisor=1.0;
            end
            w=xx([border,1],'corner')/borderdivisor;  % make a linear weigth from left to right
            if (method == 3)
                w=(1-cos(w*pi/2));  % cos ^2 weighting
            else
                w=(1-cos(w*pi))/2;  % cos ^2 weighting
            end
            wi=w(end:-1:0);
            Avg=gaussf(im);
            Avg=(Avg(0)+Avg(end))/2;
            tmp=im2(0:border-1)+Avg.*wi; 
            im2(0:border-1)=tmp;
            tmp=sumw(0:border-1)+wi;
            sumw(0:border-1)=tmp;
            tmp=max(maxw(0:border-1),wi(:));
            maxw(0:border-1)=tmp;
            tmp=im2(end-(border-1):end)+Avg.*w;
            im2(end-(border-1):end)=tmp;
            tmp=sumw(end-(border-1):end)+w;
            sumw(end-(border-1):end)=tmp;
            tmp=max(maxw(end-(border-1):end),w(:));
            maxw(end-(border-1):end)=tmp;
        else
            border=ceil(size(im,1)/2*(part/2));
            if (border > 1)
                if (method>=3)
                    borderdivisor=(border-0.5);  % This is used for the right scaling of the weights near the edges
                else
                    borderdivisor=(border-1);  % For compatibility reasons this is kept (method 2). This means that the left and right pixels are both identical
                end
            else
                borderdivisor=1.0;
            end
            w=xx([border,size(im,2)],'corner')/borderdivisor;  % make a linear weigth from left to right
            if (method == 3)
                w=(1-cos(w*pi/2));  % cos ^2 weighting
            else
                w=(1-cos(w*pi))/2;  % cos ^2 weighting
            end
            wi=w(end:-1:0,:);
            % Avg=repmat(permute(expanddim(gaussf(squeeze(im(0,:)/2+im(end,:)/2),sigmaY),2),[2 1]),[border 1]);
            Avg=permute(expanddim(gaussf(squeeze(im(0,:)/2+im(end,:)/2),sigmaY),2),[2 1]);
            tmp=im2(0:border-1,:)+Avg.*wi; 
            im2(0:border-1,:)=tmp;
            tmp=sumw(0:border-1,:)+wi;
            sumw(0:border-1,:)=tmp;
            tmp=max(maxw(0:border-1,:),wi);
            maxw(0:border-1,:)=tmp;
            tmp=im2(end-(border-1):end,:)+Avg.*w;
            im2(end-(border-1):end,:)=tmp;
            tmp=sumw(end-(border-1):end,:)+w;
            sumw(end-(border-1):end,:)=tmp;
            tmp=max(maxw(end-(border-1):end,:),w);
            maxw(end-(border-1):end,:)=tmp;
        end
        if (dim >= 2)
            border=ceil(size(im,2)/2*(part/2));
            if (border > 1)
                if (method>=3)
                    borderdivisor=(border-0.5);  % This is used for the right scaling of the weights near the edges
                else
                    borderdivisor=(border-1);  % For compatibility reasons this is kept (method 2). This means that the left and right pixels are both identical
                end
            else
                borderdivisor=1.0;
            end
            w=yy([size(im,1),border],'corner')/borderdivisor;  % make a linear weigth from left to right
            if (method == 3)
                w=(1-cos(w*pi/2));  % cos ^2 weighting
            else
                w=(1-cos(w*pi))/2;  % cos ^2 weighting
            end
            wi=w(:,end:-1:0);
            Avg=gaussf(squeeze(im(:,0)/2+im(:,end)/2),sigmaX); % repmat(gaussf(squeeze(im(:,0)/2+im(:,end)/2),sigmaX),[1 border])
            tmp=im2(:,0:border-1)+Avg.*wi; 
            im2(:,0:border-1)=tmp;
            tmp=sumw(:,0:border-1)+wi;
            sumw(:,0:border-1)=tmp;
            tmp=max(maxw(:,0:border-1),wi);
            maxw(:,0:border-1)=tmp;
            tmp=im2(:,end-(border-1):end)+Avg.*w;
            im2(:,end-(border-1):end)=tmp;
            tmp=sumw(:,end-(border-1):end)+w;
            sumw(:,end-(border-1):end)=tmp;
            tmp=max(maxw(:,end-(border-1):end),w);
            maxw(:,end-(border-1):end)=tmp;
        end
        tmp=im2(sumw>0)./sumw(sumw>0);
        im2(sumw>0)=tmp;
        de=im2*maxw+im*(1-maxw);
    elseif ndims(im) == 4
        de=im*0;noiseFactor=im*0;  % The zeros lead to a copy and avoid a problem in cudaMat
        for e=0:size(im,4)-1
            myelem=squeeze(im(:,:,:,e));
            [de2,noiseFactor2]=DampEdge (myelem,part,dim,shape,method,sigma);
            de(:,:,:,e)=de2;
            noiseFactor(:,:,:,e)=noiseFactor2;
        end
        return
    else
        error('DampEdges: Maximally four-dimensional images supported for mirror-eges')
    end
else
        error('DampEdges: Unknown method. Only 1 and 2 are allowed.')    
end
dipsetpref('DerivativeFlavour',mypref)
de=expanddim(de,origdims);

if nargout > 1
    noiseFactor = sum((1-maxw).^2.*im)/sum(im);  % estimates the reduction in noise power when measured in Fourier space
    % This equation takes the square since the influence on the noise is
    % according to the variance of I, which is proportional to the square
    % of I
end
if nargout > 2
    offsetFactor = sum((1-maxw).^2.)/prod(size(maxw));  % estimates the reduction in noise power when measured in Fourier space
end

% if length(sizevector3D)>2 & dim==2
%     i1=newim(sizevector3D);
%     i1(:,:,0)=im(:,:);
%     im=i1;
% end
