% out=extract(img,asize,center,value,blurSize, blurStart) : will extract (cut, pad) a specified region centered at a defined. Also with padding/extrapolation/inpainting.
% img: input image
% asize: size of output image (default:[] = original size)
% center: center position (will be rounded) (default : [] = center of image)
% value: value to use for padding (default: 0)
%        if 'cyclic' is specified for value a cyclic blending ist done
%        if 'gcyclic' is specified for value a cyclic blending ist done with Gaussian prefiltering
% blurSize: The maximal kernelsize to be used for 'gcyclic'. default: 8.0;
% blurStart: Extra number of pixels to start blurring earlier (similar to DampEdge). default: [] = zeros;
%
% Example:
% a  = readim('orka');
% q = extract(a + 0.0,[300 300],[],'gcyclic')
% cat(1,ft(a),extract(ft(q),[256 256]))
%
% Author: Rainer Heintzmann
% Revision: 05/02/2017: Not using "eval" any more and Gaussian behaviour was improved.
%
% See also: DampEdge
function out=extract(img,asize,center,value,blurSize,blurStart)
global use_newim_cuda
isize=size(img);
if nargin < 4
    value=0;
end
if nargin < 3 || isempty(center)
    center=floor(isize/2);
else
end
idim=size(asize);
if idim(1) > 1
        asize = asize';
end
idim=size(center);
if idim(1) > 1
        center = center';
end
if length(asize) < length(isize)
    for d=length(asize)+1: length(isize)
        asize(d) = size(img,d);
    end
end
if length(center) < length(isize)
    for d=length(center)+1: length(isize)
        center(d) = floor(size(img,d)/2);
    end
end
N=ndims(img);
if isa(img,'cuda')
    out=newim(asize,datatype(img));   % problems with type conversion: +value;
else
    tmp=use_newim_cuda;
    use_newim_cuda=0;
    out=newim(asize,datatype(img));   % problems with type conversion: +value;
    use_newim_cuda=tmp;    
end

srccenter=round(center');
if length(asize) > length(srccenter)
   srccenter(end+1:length(asize))=0;
   isize(end+1:length(asize))=1;
end

srcstart=srccenter-floor(asize'/2);
srcend=srcstart+asize'-1;
dststart=zeros(length(asize),1);
dststart(srcstart<0)=-srcstart(srcstart<0);
dstend=asize-1;
dstend(srcend>=isize')=dstend(srcend'>=isize)'-srcend(srcend'>=isize)+isize(srcend'>=isize)'-1;
srcend(srcend>=isize')=isize(srcend'>=isize)-1;
srcstart(srcstart<0)=0;

if isnumeric(value)  
    if value~=0
        out=out+value;  % the rest will be overwritten below
    end
end

% s = 'out(';  % out version using "eval". Now replaced by subsasgn and subsref
% for ii=1:ndims(out)
%    s = [s 'dststart(' num2str(ii) '):dstend(' num2str(ii) '),'];
% end
% s(end) =[];s = [s ')=img('];
% for ii=1:ndims(img)
%    s = [s 'srcstart(' num2str(ii) '):srcend(' num2str(ii) '),'];
% end
% s(end) =[];s = [s ');'];
% eval(s); 

idxS=num2cell(repmat(':',[1 ndims(img)]));idxD=num2cell(repmat(':',[1 ndims(img)]));
for ii=1:ndims(img)
    if ii<=numel(srcstart)
        idxS{ii}=[srcstart(ii):srcend(ii)];
    end
    if ii<=numel(dststart)
        idxD{ii}=[dststart(ii):dstend(ii)];
    end
end
S=struct('type','()','subs',{idxS}); D=struct('type','()','subs',{idxD});
out=subsasgn(out,D,subsref(img,S));   % performs the assignment

if ~isnumeric(value)  % other case was dealt with above
    switch value
        case {'cyclic','gcyclic'}
            for d=1:ndims(out) % dimension
                pstart=dststart(d);
                pend=dstend(d);
                sges=pstart+(size(out,d)-pend);  % total size to fill
                p1 = SubSlice(out,d,pstart);
                p2 = SubSlice(out,d,pend);
%                 eval(['p1=out' hyperplaneString(d,pstart,ndims(out)) ';']);
%                 eval(['p2=out' hyperplaneString(d,pend,ndims(out)) ';']);
                for p=0:pstart-1 % position
                    w=(p+(size(out,d)-pend))/sges;
                    w=(sin(w*pi/2))^2;
                    aline=p1*w+p2*(1-w);
                    %  eval(['out' hyperplaneString(d,p,ndims(out)) '=aline;']);
                    idxD=num2cell(repmat(':',[1 ndims(out)]));
                    idxD{d}=p; D=struct('type','()','subs',{idxD});
                    out=subsasgn(out,D,aline);   % performs the assignment                    
                end
                for p=pend+1:size(out,d)-1 % position
                    w=(p-pend)/sges;
                    w=(sin(w*pi/2))^2;
                    aline=p1*w+p2*(1-w);
                    % eval(['out' hyperplaneString(d,p,ndims(out)) '=aline;']);
                    idxD=num2cell(repmat(':',[1 ndims(out)]));
                    idxD{d}=p; D=struct('type','()','subs',{idxD});
                    out=subsasgn(out,D,aline);   % performs the assignment
                end
            end
            if value(1)=='g'  % use gaussian blurring: run through the regions again in a 2nd pass and do the blurring
                if nargin < 5
                    blurSize = 8.0;
                end
                if nargin < 6
                    blurStart = zeros(1,ndims(out));
                else
                    if numel(blurStart)==1
                        blurStart = repmat(blurStart,[1 ndims(out)]);
                    end
                end
               for d=1:ndims(out) % dimension
                pstart=dststart(d)+blurStart(d);
                pend=dstend(d)-blurStart(d);
                sges=pstart+(size(out,d)-pend);  % total size to fill
                p1 = SubSlice(out,d,dststart(d));
                p2 = SubSlice(out,d,dstend(d));
                for p=0:pstart-1 % position
                    w=(p+(size(out,d)-pend))/sges;
                    w=(sin(w*pi/2))^2;
                    aline=p1*w+p2*(1-w);
                    aline=gaussf(aline,w*(1-w)*blurSize);  % The amount of blurring varies 
                    idxD=num2cell(repmat(':',[1 ndims(out)]));
                    idxD{d}=p; D=struct('type','()','subs',{idxD});
                    p3 = SubSlice(out,d,p);
                    w2 = abs(w-0.5) * 2.0;
                    out=subsasgn(out,D,aline .*(1.0-w2) + p3 .* w2);   % performs the assignment                    
                end
                for p=pend+1:size(out,d)-1 % position
                    w=(p-pend)/sges;
                    w=(sin(w*pi/2))^2;
                    aline=p1*w+p2*(1-w);
                    aline=gaussf(aline,w*(1-w)*blurSize); % The amount of blurring varies 
                    % eval(['out' hyperplaneString(d,p,ndims(out)) '=aline;']);
                    idxD=num2cell(repmat(':',[1 ndims(out)]));
                    idxD{d}=p; D=struct('type','()','subs',{idxD});
                    p3 = SubSlice(out,d,p);
                    w2 = abs(w-0.5) * 2.0;
                    out=subsasgn(out,D,aline .*(1.0-w2) + p3 .* w2);   % performs the assignment
                end
               end
            end
        otherwise
            error('extract: Unknown blend method');
    end
        
end
    