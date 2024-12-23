global KhorosOutput %KhorosOuput is returned to this global variable
global KhorosRoot

mainpath = 'C:\Users\Dina Ratsimandresy\Documents\GitHub\PSFToolbox-Matlab\'; % CHANGE DIRECTORY HERE
cd(mainpath);
addpath(genpath([mainpath,'\Code\']));

KhorosRoot=[mainpath '\Code\khorosBin\'];
global KhorosOutput %KhorosOuput is returned to this global variable

DIPPath='C:\Program Files\DIPimage 2.9'; 
addpath([DIPPath '\common\dipimage']);
dip_initialise;
dipsetpref('ImageFilePath', [DIPPath '\images']);
dipsetpref('DefaultMappingMode','lin')
dipsetpref('TrueSize','off')
set(0,'DefaultFigurePaperType','A4');

%% Parameters
ImageParam=struct('Sampling',[65.5 65.5 160],'Size',[512 512 32]);
PSFParam=struct('NA',1.4,'n',1.518,'lambdaEm',542,'Aplanatic',-1);
AddParams=struct('ns',1.33,'ng',1.518,'ni',1.515,'ng0',1.518,'ni0',1.518,'ts',1e3,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5,'tilt',[0.6 0]);
MethodList={"RichardsWolff", 'RichardsWolffInt','SlicePropagation',"CZT", 'VolumeShell', 'SincR'};
%% Compute PSFs
hAll=[];% store PSFs
for k=1:1:length(MethodList)
    tic; h=GenericPSFSim(ImageParam,PSFParam,MethodList{k},AddParams); t1=toc; % aberrated
    hAll=cat(4,hAll,h);
end

hAll % display