function [R, wls, zpos, ts, datainfo] = loadMSOTRecon(varargin)
% LOADMSOTRECON  Load MSOT Reconstructions
%   [R wls zpos ts datainfo] = LOADMSOTRECON()       
%           opens Open File Dialog to choose file, then loads the first
%           Recon
%   [R wls zpos ts datainfo] = LOADMSOTRECON(file)  
%           opens the specified file, then loads the first recon
%   [R wls zpos ts datainfo] = LOADMSOTRECON(datainfo)   
%           opens the specified file, then loads the first recon
%   [R wls zpos ts datainfo] = LOADMSOTRECON(datainfo,reconNum)   
%           opens the specified file, then loads the specified recon
%   [R wls zpos ts datainfo] = LOADMSOTRECON(datainfo,reconNum,selMat)   
%           opens the specified file, then loads the specified recon
%           restricted to selMat, which is a multidimensional structure of 
%           frame ids. For all, use datainfo.ReconNode(i).ReconStructure
% 
% This function loads the MSOT Reconstructions from the selected scan
% folder. 
%
%
% Return values:
% R          An array containing the reconstructed images in 8D
%            1) x
%            2) y
%            3) RUN (repetition of the whole acquisition step, time
%               dimension (timepoints as separate vector)
%            4) z-position (if 2D system with translation stage)
%            5) REPETITION (currently unused)
%            6) Wavelength (see wls vector)
%            7) Individual Images (only if not averaged)
% wls        Wavelength Vector
% zpos       Array of z-stage positions
% ts         Multidimensional array of timestamps in s, same 
%            dimension as R

%% parameters
rNode = 0;
selMat = [];
datainfo = [];

% if no arguments given
if (nargin == 0 || isempty(varargin{1}))
    datainfo = loadMSOT();
else
    datainfo = varargin{1};
end
if ~isstruct(datainfo) && ischar(datainfo)
    datainfo = loadMSOT(datainfo);
end

% if no reconNode numner is given, use the first
if nargin >= 2
    rNode = varargin{2};
else
    rNode = 1;
end

if (numel(datainfo.ReconNode) < rNode)
    error('Reconstruction not found');
end

% if no selMat is given, use all
if nargin >= 3 && ~isempty(varargin{3})
    selMat = varargin{3};
else
    selMat = datainfo.ReconNode(rNode).ReconStructure;
end

%% Initialise Progress Bar
wbar = waitbar(0,'Initialising, please wait...' ) ;


%% initialise return array
sel = ~isnan(selMat);   % this will retain structure, marking only not-nans
svec = size(sel);       % retain structural information
selMat = selMat(sel);   % remove nans and make selMat 1D

% image return vector
n = datainfo.ReconNode(rNode).Resolution;
R = nan([n n prod(svec)],'double');
% initial vectors for wavelengths and zpositions and timestamps
wls = zeros(size(selMat));
zpos = zeros(size(selMat));
ts = zeros(prod(svec),1);

%% Load Images
fname = [datainfo.RealPath '\RECONs\' datainfo.ReconNode(rNode).GUID '.bin'];
[FID str] = fopen(fname,'r','n');
if (FID == -1)
    error(['Error opening binary Recon File ' fname ': ' str]);
    return;
end
% handle progress bar
waitbar(0,wbar,'Loading Reconstructions...');

% load all selected images sequentially
for j = 1:numel(selMat)
    waitbar(j/numel(selMat),wbar);
    id = double(selMat(j));
    if id == 0 || isnan(id)
        fprintf('Skipping frame %i - id is 0\n',j);
        continue;
    end
    
    % copy Meta-Information from Scan Frame
    wls(j) = datainfo.ReconNode(rNode).Frames(id).Wavelength;
    zpos(j) = datainfo.ReconNode(rNode).Frames(id).ZPos;
    ts(j) = datainfo.ReconNode(rNode).Frames(id).RelTime;
    try
        startbyte = (id-1)*n*n*8;    %determine start position (double prec)
        fseek(FID,startbyte,-1);     % move towards position
        R(:,:,j) = fread(FID,[n n],'double');
    catch ex
        warning(['Cannot Read ReconFrame ' num2str(id) ', skipping...']);
    end
   
end
fclose(FID);

%% reshape data to fit requirements
R = reshape(R,[n n svec]);
ts = reshape(ts,svec);
wls = unique(wls);
zpos = unique(zpos);

close(wbar) ;
