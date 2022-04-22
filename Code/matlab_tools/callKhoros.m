%callKhoros  calls a Khoros routine and returns the output images
%
% callKhoros(fktname, in, arglist)
% 
% EXAMPLE:
% callKhoros('kimport',{i1,i2},{'wsize',256;'hsize',256;'filename','bla'})

% (C) Copyright 2004           Department of Molecular Biology
%     All rights reserved      Max-Planck-Institute for Biophysical Chemistry
%                              Am Fassberg 11, 37077 G"ottingen
%                              Germany
%
% Bernd Rieger & Rainer Heintzmann, June, 2004.


%global KhorosOutput from startup file
function varargout = callKhoros(fktname, in, arglist)
if nargin ~=3
   error('Usage: [out1,..] = callKhoros(fktname, in, arglist).');
end
if isempty(tempdir)
   error('OS does not give a temp directory.');
end
if  ~iscell(arglist) | ~iscell(in)
   error('inputs not cell arrays');
end
if size(arglist,2)~=2
   error('arglist must be of form {''tag1'',256;''tag2'',256}')
end

td = tempdir;
kh = [fktname ' '];

[s,sysid]=system('echo $$');
sysid(end)=[];%remove \n
nlist = size(arglist,1);

%Remove '__input' and '__output'
for ii=1:nlist
   if isstr(arglist{ii,2})
      if strcmp(arglist{ii,2},'__ouput') | strcmp(arglist{ii,2},'__input')
         arglist{ii,:}=[]; %remove it from this from the list
      end
   end
end
nlist = size(arglist,1);
%nlist

for ii=1:length(in)
   if isstr(in{ii})
     fn_in{ii} = in{ii};
   else
     %fn_in{ii}= [td 'xxx']; 
     fn_in{ii}= [td  'MatlabKhoros_in_' sysid '.' num2str(ii)];
   end
end

%Add parameters 
no=0;
ni=0;
for ii=1:nlist
if ~strcmp(arglist{ii,1},'')
   if isstr(arglist{ii,2})
      if strcmp(arglist{ii,2},'OK_out')
         no = no+1;
         fn_out{no}= [td  'MatlabKhoros_out_' sysid '.' num2str(no)];
         %fn_out{no}= [td  'xxx'];
         kh = [kh '-' arglist{ii,1} ' ' fn_out{no} ' '];
         
      elseif strcmp(arglist{ii,2},'OK_in')
         ni = ni+1;
         kh = [kh '-' arglist{ii,1} ' ' fn_in{ni} ' '];
         
      else
         kh =[kh '-' arglist{ii,1} ' ' arglist{ii,2} ' '];
      end
   else
      kh =[kh '-' arglist{ii,1} ' ' num2str(arglist{ii,2}) ' '];
   end
end
end

%Write input files to disk
%in
%fn_in

for ii=1:ni
   if ~isstr(in{ii})
       if isempty(in{ii})
           error(sprintf('Trying to pass an empty image to khoros as argument number %d.',ii));
       end
     writekhoros(in{ii},fn_in{ii});
   end
end

%call Khoros
fprintf(' Khoros called with command:\n %s\n',kh)
global KhorosOutput
if strcmp(KhorosOutput,'Direct')
    s=system(kh);
else    
    [s,KhorosOutput]=system(kh);
end

if s
   error(['Khoros command failed with: ' KhorosOutput]);
end

%Read output images from disk,delete and return them
if nargout < no
    fprintf('WARNING: %s: Only %d out of %d output arguments assigned.\n',fktname,nargout,no);
    no=nargout;
end
if nargout > no
    error('Too many output arguments\n');
end

for ii=1:no
   %varargout{ii} =1;
   varargout{ii}=readkhoros(fn_out{ii});
   if ispc
   system(['del ' fn_out{ii}]);
   else
   system(['/bin/rm ' fn_out{ii}]);
   end
end

for ii=1:ni
  if ~isstr(in{ii})
   if ispc
   system(['del ' fn_in{ii}]);
   else
   system(['/bin/rm ' fn_in{ii}]);
   end
  end
end
