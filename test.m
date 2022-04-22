
% mainpath = 'C:\Users\labuser\OneDrive - Stellenbosch University\PhDFolder\PSFToolboxAb\';
mainpath = 'C:\Users\Dina\OneDrive - Stellenbosch University\PhDFolder\PSFToolboxAb\';
cd(mainpath)
addpath(genpath([mainpath,'\Code\ElectricFields\']))
addpath(genpath([mainpath,'\Code\ElectricFields\psfModels']))
addpath([mainpath,'Code\ElectricFields\psfModels\mex'])
addpath(genpath([mainpath,'\Code\khoros\']))
addpath(genpath([mainpath,'\Code\matlab_tools\']))
%% Parameters
ImageParam=struct('Sampling',[65.5 65.5 160],'Size',[512 512 32]);
PSFParam=struct('NA',1.4,'n',1.518,'lambdaEm',542,'Aplanatic',-1);
AddParams=struct('ns',1.33,'ng',1.518,'ni',1.515,'ng0',1.518,'ni0',1.518,'ts',1e3,'tg',1.7e5,'tg0',1.7e5,'wd',1.5e5,'tilt',[0.6 0]);
MethodList={'RichardsWolffInt','SlicePropagation','VolumeShell'};
%% Compute PSFs
hAll=[];% store PSFs
for k=1:1:length(MethodList)
    tic; h=GenericPSFSim(ImageParam,PSFParam,MethodList{k},AddParams); t1=toc; % aberrated
    hAll=cat(4,hAll,h);
end

hAll % display