function [sig, datainfo] = loadMSOTSignals(varargin)
% LOADMSOTSIGNALS  Load binary scan data from .bin 
%   [sigMat, datainfo] = LOADMSOTSIGNALS()  
%        uses Open File Dialog to choose .msot file
%   [sigMat, datainfo] = LOADMSOTSIGNALS(filename)  
%        uses filename and path .msot file to load Data
%   [sigMat, datainfo] = LOADMSOTSIGNALS(datainfo)  
%        load signals from the scan described by datainfo
%   [sigMat, datainfo] = LOADMSOTSIGNALS(filename/datainfo,selMat)  
%        restricts loading of data to specific frame numbers supplied in
%        selMat. If empty, datainfo.ScanStructure is used. Output
%        will have the same structure as selMat.
%   [sigMat, datainfo] = LOADMSOTSIGNALS(filename/datainfo,selMat,params)  
%        supplies an additional parameter struct with the following fields:
%        - average: Average data while loading (default: 1)
%              NOTE: Only applies if data was acquired without averaging
%        - usePower: Laser Energy Handling (default: 1)
%            * 0: load plain binary information (will be scaled to fit
%                 16bit). IMPORTANT: Only use this for non-averaged data!
%            * 1: undo scaling and retain laser energy correction
%            * 2: undo scaling and laser energy correction.
%                 WARNING: requires subsequent laser energy correction
%        - shotsel: Array of single shots to average. Empty: all (default)
%              NOTE: Only applies on non-averaged data.
%        - filter: Filter the data after loading (default: 0)
%            * 0: Do not filter
%            * 1: Filter (zero phase filtering, i.e. forward and backward)
%            * 2: Filter (non-zero phase filtering, i.e. only forward)
%        - f_LPF: Low-pass edge frequency in Hz (default: 8e6)
%        - f_HPF: High-pass edge frequency in Hz (default: 50e3)
%        - f_HPFOrder: Order of Chebychev I HP Filter (default: 4)
%        - f_HPFRipple: Allowed passband ripple at edge freq (default: 0.1)
%
%
% Return Values:
%   1: sigMat   a matrix with dimensions equivalent to
%               datainfo.ScanStructure or selMat if supplied
%
% Example:
%   params.filter = 1;          % filter signals
%   sigmat = loadMSOTSignals('Scan_1\Scan_1.msot',[],params);
% 


selMat = [];
sig = [];
datainfo = [];
par.shotsel = [];           % select shots to load, empty is all
par.usePower = 1;        % undo scaling or apply power correction
par.f_LPF = 8000000;        % filter freq. defaults
par.f_HPF = 50000;          % filter freq. defaults
par.f_HPFOrder = 4;
par.f_HPFRipple = 0.1;
par.filter = 0;             % dont filter, 1: zerophase, 2:non-zp
par.average = 1;            % but average per default

fs = 4e7;

% get function parameters
if nargin == 0
    datainfo = loadMSOT();
end
if (nargin >= 1)
    if (isstruct(varargin{1}))
        datainfo = varargin{1};
    else
        datainfo = loadMSOT(varargin{1});
    end
end
if nargin >= 2
    selMat = varargin{2};
end
if nargin >= 3
    cpar = varargin{3};
    % transfer all fields to parameter array
    fx = fieldnames(cpar);
    for j = 1:numel(fx)
        par = setfield(par,fx{j},getfield(cpar,fx{j}));
    end
    clear cpar j fx;
end

% check parameters
if (isempty(datainfo)), 
    warning('Metadata not found, exiting...');
    return;
end
% per default, load complete dataset
if (isempty(selMat)), selMat = datainfo.ScanStructure; end

%% Progress Bar
wbar = waitbar(0,'Opening file...');

%% open File
FID = fopen(datainfo.FileName,'r');

%% determine read structure

% select shots
if ( ~isempty(par.shotsel) && ndims(selMat) >= 5 )
   selMat = selMat(:,:,:,:,par.shotsel);
end


svec = size(selMat);        % dimension vector
% if not averaged, read structure as is in 1D
if ~par.average
    frameids = reshape(selMat,prod(svec),1);
    numavg = 1;
% otherwise put ids to average in 2D array
else
    % use last (5th) number for number of averages and shorten dim vector
    if numel(svec) == 5,    
        numavg = svec(5); 
        svec = svec(1:4);       % cut size vector for output
    % otherwise the dataset is already averaged
    else numavg = 1; end;
    
    frameids = reshape(selMat,prod(svec),numavg);
end

%% data type
dtype = 'uint16';                           % default to native uint16
if numavg > 15  dtype = 'uint32'; end       % for more averages use uint32
% if filtering or power correction required, double is necessary
if (par.filter || par.usePower) dtype = 'double'; end
    
%% initialize filters
if (par.filter)
    if par.f_LPF, [b_LPF,a_LPF] = cheby1( 8, .01, 2 * par.f_LPF/fs * .9 ) ; end
    
    if par.f_HPF, 
        if par.f_HPFOrder == 4
            [b_HPF,a_HPF] = cheby1( 4, .01, 2 * par.f_HPF/fs * 1.46, 'high' ) ; 
        elseif par.f_HPFOrder == 1
            [b_HPF,a_HPF] = cheby1( 1, 1, 2 * par.f_HPF/fs , 'high' ) ; 
        end
    end
end

%% read data
framenum = numel(frameids)-1;
waitbar(0,wbar,'Loading Signals...');
numproj = datainfo.HWDesc.NumDetectors;
numsamples = datainfo.MeasurementDesc.RecordLength;
sig = zeros(numsamples,numproj,size(frameids,1),dtype);

for j = 1:size(frameids,1);
    waitbar(j/framenum,wbar);
    tmp = zeros(numsamples,numproj,dtype);
    readcount = 0;          % counter for read frames
    % iterate through second dimension to capture averages
    for ii = 1:size(frameids,2)
        % get current id and frame data
        id = frameids(j,ii);
        if (isnan(id)) continue; end;
        frame = datainfo.ScanFrames(id);
        offset = frame.IDOffset;
%         fprintf('id: %i - offset: %i\n',id,offset);
        
        readcount = readcount + 1;

        % position and read uint16 data
        tmp2 = zeros(numsamples,numproj,dtype);
        try
            fseek(FID,(offset)*numsamples*numproj*2,-1);
            tmp2(:,:) = fread(FID,[numsamples numproj],'uint16');
        catch ex
            warning(['Cannot Read ScanFrame ' num2str(id) ', skipping...']);
            readcount = readcount - 1;
        end
        
        % correct for scaling
        if par.usePower
            if (par.usePower == 2)
                tmp2 = tmp2 ./ frame.CorrectionFactor .* frame.LaserEnergy;
            else
                tmp2 = tmp2 ./ frame.CorrectionFactor;
            end
        end
        
        % filtering
        if par.filter == 1
            if par.f_LPF, tmp2 = FiltFiltM( b_LPF, a_LPF, tmp2, 1, 2 ) ; end
            if par.f_HPF, tmp2 = FiltFiltM( b_HPF, a_HPF, tmp2, 1, 2 ) ; end
        elseif par.filter == 2
            if par.f_LPF, tmp2 = FilterM( b_LPF, a_LPF, tmp2 ) ; end
            if par.f_HPF, tmp2 = FilterM( b_HPF, a_HPF, tmp2 ) ; end
        end
        
        % add to vector
        tmp = tmp + tmp2;
        
    end
    % average all if any frames were read (otherwise leave zeros)
    if readcount
        sig(:,:,j) = tmp ./ readcount;
    end
end;
sig = reshape(sig,[numsamples numproj svec]);
clear ii j tmp tmp2 frame id;

% hide dialog
close(wbar);

%% close file
fclose(FID);