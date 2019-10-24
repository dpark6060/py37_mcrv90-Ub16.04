function [R, zpos, ts, datainfo] = loadMSOTMsp(varargin)
% LOADMSOTRECON  Load MSOT Reconstructions
%   [R wls zpos ts datainfo] = LOADMSOTMSP()       
%           opens Open File Dialog to choose file, then loads the first
%           Msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(file)  
%           opens the specified file, then loads the first msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(datainfo)   
%           opens the specified file, then loads the first msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(datainfo,mspNum)   
%           opens the specified file, then loads the specified msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(datainfo,mspNum,selMat)   
%           opens the specified file, then loads the specified msp
%           restricted to selMat, which is a multidimensional structure of 
%           frame ids. For all, use datainfo.MspNode(i).Structure
% 
% This function loads the MSOT MSPs from the selected scan folder. 
%
%
% Return values:
% R          An array containing the MSPed images in 8D
%            1) x
%            2) y
%            3) z (only if 3D detector)
%            4) RUN (repetition of the whole acquisition step, time
%               dimension (timepoints as separate vector)
%            5) z-position (if 2D system with translation stage)
%            6) REPETITION (currently unused)
%            7) Wavelength (see wls vector)
%            8) Individual Images (only if not averaged)
% wls        Wavelength Vector
% zpos       Array of z-stage positions
% ts         Multidimensional array of timestamps in s, same 
%            dimension as R

%% parameters
mNode = 0;
selMat = [];
datainfo = [];

% if no arguments given
if (nargin == 0 || isempty(varargin{1}))
    datainfo = loadMSOT();
else
    datainfo = varargin{1};
end
if ~isstruct(datainfo) || isempty(datainfo)
    return
end

% if no reconNode numner is given, use the first
if nargin >= 2
    mNode = varargin{2};
else
    mNode = 1;
end

if (numel(datainfo.MSPNode) < mNode)
    error('MSP not found');
end

% if no selMat is given, use all
if nargin >= 3
    selMat = varargin{3};
else
    selMat = datainfo.MSPNode(mNode).Structure;
end

%% Initialise Progress Bar
wbar = waitbar(0,'Initialising, please wait...' ) ;


%% initialise return array
sel = ~isnan(selMat);   % this will retain structure, marking only not-nans
svec = size(sel);       % retain structural information
selMat = selMat(sel);   % remove nans and make selMat 1D

% image return vector
n = datainfo.MSPNode(mNode).Resolution;
R = nan([n n prod(svec)],'double');
% initial vectors for wavelengths and zpositions and timestamps
wls = zeros(size(selMat));
zpos = zeros(size(selMat));
ts = zeros(prod(svec),1);

%% Load Images
fname = [datainfo.RealPath '\MSPs\' datainfo.MSPNode(mNode).GUID '.bin'];
[FID str] = fopen(fname,'r','n');
if (FID == -1)
    error(['Error opening binary MSP File ' fname ': ' str]);
    return;
end
% handle progress bar
waitbar(0,wbar,'Loading MSP...');

% load all selected images sequentially
for j = 1:numel(selMat)
    waitbar(j/numel(selMat),wbar);
    id = double(selMat(j));
    if id == 0, warning('Skipping frame...'); continue; end;
    slid = datainfo.MSPNode(mNode).SliceIndex(id);
    
    % copy Meta-Information from Scan Frame
    zpos(j) = datainfo.MSPNode(mNode).Slices(slid).ZPos;
    ts(j) = datainfo.MSPNode(mNode).Slices(slid).RelTime;
    try
        startbyte = (id-1)*n*n*8;      % determine start position (double prec)
        fseek(FID,startbyte,-1);     % move towards position
        R(:,:,j) = fread(FID,[n n],'double');
    catch ex
        warning(['Cannot Read MspFrame ' num2str(id) ', skipping...']);
    end
   
end
fclose(FID);

%% reshape data to fit requirements
R = reshape(R,[n n svec]);
ts = reshape(ts,svec);
zpos = unique(zpos);

close(wbar) ;
