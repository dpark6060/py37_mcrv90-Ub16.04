function datainfo  = loadMSOT( varargin )
% LOADMSOT  Load .msot META information
%   datainfo = LOADMSOT()       opens Open File Dialog to choose file
%   datainfo = LOADMSOT(file)   opens the specified file
% 
% This function loads the META information on the selected scan subfolder
% including information about both RECONs and MSPs performed on the scan
% raw data.
%
% Required Parameters:
%   <none>
%
% Optional:
%   1: filename    path and name of the file to load.
%
% Return Values:
%   1: datainfo    a structure containing the meta information, [] on error 
%
% Example:
%   datainfo = loadMSOT('Scan_1/Scan_1.msot');


%% parameters
if ( nargin == 0 )
    
    path = getappdata( 0, 'last_path' ) ;
    if( isempty( path ) )
        path = cd ;
    end ;
    fc = javax.swing.JFileChooser( java.io.File( path ) ) ;
    fc.setAcceptAllFileFilterUsed(0) ;
    fc.addChoosableFileFilter( javax.swing.filechooser.FileNameExtensionFilter(  'MSOT Header',  'msot' ) ) ;
    format = 'xml' ;
    
    if( fc.showOpenDialog( [] ) == javax.swing.JFileChooser.APPROVE_OPTION )
        path = char( fc.getSelectedFile().getParent() ) ;
        setappdata( 0, 'last_path', path ) ;
        filename = char( fc.getSelectedFile().getPath() ) ;
    else
        sigMat = [] ;
        datainfo = [] ;
        return ;
    end ;
elseif( nargin == 1 )
    filename = varargin{1} ;
else
    datainfo = [];
    error( 'Not enough input arguments' ) ;
end ;

%% Progress Bar
wbar = waitbar(0,'Initialising...');

%% FileName
fstr = java.lang.String( filename ) ;
xml_filename = '' ;
if( fstr.endsWith( '.msot' ) )
    xml_filename = filename ;
    filename = char( fstr.substring( 0, fstr.lastIndexOf( '.msot' ) ).concat( '.bin' ) ) ;
end ;
lastslash = fstr.lastIndexOf( '/' );
if lastslash < 0,
    dirname = './';
else
    dirname = char( fstr.substring( 0, lastslash ) ) ;
end

%% XML header format
waitbar(0,wbar,'Parsing XML...');
try
    patdoc = javaMethod( 'parse', 'com.itheramedical.msotbeans.DataModelMsotProjectDocument$Factory', java.io.File( xml_filename ) ) ;
catch ex
    close(wbar)
    error( ['Selected XML file could not be loaded. Check well-formedness.\n' ex.message]) ;
end ;
dm = patdoc.getDataModelMsotProject() ;


%% declaration
waitbar(0,wbar,'Reading general information...');
datainfo = struct() ;
datainfo.FileName = filename;
datainfo.XMLFileName = xml_filename;
                
%% admin data
waitbar(0,wbar,'Reading admin data...');
datainfo.FriendlyName = char(dm.getFriendlyName);
datainfo.FolderName = char(dm.getFolderName);
datainfo.RealPath = dirname;
datainfo.CreationTimeTxt = char(dm.getCreationTime.toString);
datainfo.CreationTime = datenum(datainfo.CreationTimeTxt,'yyyy-mm-ddTHH:MM:SS.FFF');
datainfo.Scientist = char(dm.getScientist);

%% measurement description
waitbar(0,wbar,'Reading measurement description...');
md = dm.getMEASUREMENTDESC;

datainfo.MeasurementDesc.Averages = md.getAVERAGESPERPROJECTION; % number of shots averaged into one frame
datainfo.MeasurementDesc.NumShots = md.getNUMBEROFFRAMES;       % number of shots per wavelength
datainfo.MeasurementDesc.NumFrames = datainfo.MeasurementDesc.NumShots / datainfo.MeasurementDesc.Averages; % Number of frames stored in raw data
datainfo.MeasurementDesc.RecordLength = md.getRECORDEDLENGTH;
datainfo.MeasurementDesc.RepRate = md.getREPETITIONRATE.getDoubleValue;
% wavelengths
datainfo.MeasurementDesc.Wavelengths = zeros(length(md.getWAVELENGTHS.getWAVELENGTHArray),1);
for i = 1:length(md.getWAVELENGTHS.getWAVELENGTHArray)
    datainfo.MeasurementDesc.Wavelengths(i) = md.getWAVELENGTHS.getWAVELENGTHArray(i-1).getIntValue;
end
datainfo.MeasurementDesc.Sequence = char(md.getSEQUENCE);

% *** Water Absorption Coefficients
try
    nwa = length(md.getWATERABSORPTIONCOEFFICIENTS.getWATERABSORPTIONCOEFFICIENTArray);
catch ex
    nwa = 0;
end
datainfo.MeasurementDesc.WaterAbsorptionCoeff = zeros(nwa,1);
for i = 1:nwa
    datainfo.MeasurementDesc.WaterAbsorptionCoeff(i) = double(md.getWATERABSORPTIONCOEFFICIENTS.getWATERABSORPTIONCOEFFICIENTArray(i-1).getCoefficient);
end
% Path Length in Water (for Water correction)
datainfo.MeasurementDesc.PathLengthInWater = nan;
try
    datainfo.MeasurementDesc.PathLengthInWater = md.getPATHLENGTHINWATER.getDoubleValue;
end

% *** Average Laser Energy Table (used if Average Energy Correction applied)
try
    nfle = length(md.getAVERAGEENERGYTABLE.getAVRAGEENERGYArray);
catch ex
    nfle = 0;
end
datainfo.MeasurementDesc.FactoryLaserEnergyTable = zeros(nfle,1);
for i = 1:nfle
    datainfo.MeasurementDesc.FactoryLaserEnergyTable(i) = md.getAVERAGEENERGYTABLE.getAVRAGEENERGYArray(i-1).getDoubleValue;
end

% Type of Laser Energy Correction
datainfo.MeasurementDesc.EnergyNormalisation = char(md.getLEnergyNormalization);

%% hardware description
waitbar(0,wbar,'Reading hardware description...');
hw = dm.getHARDWAREDESC;

datainfo.HWDesc.SamplingFrequency =  hw.getSAMPLINGFREQUENCY.getDoubleValue*1e6;
datainfo.HWDesc.TransducerType = char(hw.getTRANSDUCER);
datainfo.is3D = strcmp(datainfo.HWDesc.TransducerType,'msot3');
datainfo.is2D = strcmp(datainfo.HWDesc.TransducerType,'msot2');

% 2D arrays
eq = hw.getFRAMEDESC.getEQUAL;
if (~isempty(eq))
    datainfo.HWDesc.StartAngle = double(eq.getSTART);
    datainfo.HWDesc.StepAngle = double(eq.getSTEP);
    datainfo.HWDesc.EndAngle = double(eq.getEND);
    datainfo.HWDesc.NumDetectors = round((datainfo.HWDesc.EndAngle - datainfo.HWDesc.StartAngle) / datainfo.HWDesc.StepAngle)+1;
    for j = 0:numel(eq.getCONSTANTArray)-1
       cst = eq.getCONSTANTArray(j);
       if cst.getAxisRef == 2
           datainfo.HWDesc.Radius = cst.getDoubleValue;
       elseif cst.getAxisRef == 3
           datainfo.HWDesc.RadiusZ = cst.getDoubleValue;
       end
    end
else
    datainfo.HWDesc.StartAngle = [];
    datainfo.HWDesc.StepAngle = [];
    datainfo.HWDesc.EndAngle = [];
    datainfo.HWDesc.Radius = [];
    datainfo.HWDesc.RadiusZ = [];
end


%3D cup
pa = hw.getFRAMEDESC.getPROJECTIONArray;
if (~isempty(pa))
   datainfo.HWDesc.NumDetectors = numel(pa);
   for j = 1:numel(pa)
       va = pa(j).getVALUEArray;
       % 3 coordinates
       datainfo.HWDesc.ProjectionData(j,:) = [va(1).doubleValue va(2).doubleValue va(3).doubleValue ];
       datainfo.HWDesc.Radius = mean(sqrt(sum(datainfo.HWDesc.ProjectionData.^2,2)));
   end
end
    
% prj = hw.getFRAMEDESC.get


%% Scan Node information
waitbar(0,wbar,'Reading scan description...');
sn = dm.getScanNode;
% wavelengths
datainfo.Wavelengths = zeros(length(sn.getWavelengths.getWavelengthArray),1);
for i = 1:length(sn.getWavelengths.getWavelengthArray)
    datainfo.Wavelengths(i) = sn.getWavelengths.getWavelengthArray(i-1).getIntValue;
end
% name and comment
datainfo.Name = char(sn.getName);
datainfo.Comment = char(sn.getComment);



%% Scan Frame Information
waitbar(0,wbar,'Reading scan frames...' ) ;
fa = sn.getScanFrames.getDataModelScanFrameArray;
maxRep = 0;
maxRun = 0;
zposArr = [];
tempavg = 0;
datainfo.ScanFrames(length(fa)) = struct();
for i = 1:length(fa)
    waitbar(i/length(fa),wbar,'Reading scan frames...' ) ;
    fr = fa(i).getFrame; 
    datainfo.ScanFrames(i).IDOffset = fa(i).getIDOffset;
    datainfo.ScanFrames(i).Temperature = fr.getTEMPERATURE.getDoubleValue;
    tempavg = tempavg+datainfo.ScanFrames(i).Temperature;
    datainfo.ScanFrames(i).Wavelength = fr.getWAVELENGTH.getIntValue;
    datainfo.ScanFrames(i).CorrectionFactor = fr.getPOWER.getDoubleValue;
    datainfo.ScanFrames(i).ZPos =  fr.getZPOS.getDoubleValue;
    datainfo.ScanFrames(i).Timestamp = fr.getTimestamp.doubleValue * 1e-7;
    datainfo.ScanFrames(i).RelTime = datainfo.ScanFrames(i).Timestamp - datainfo.ScanFrames(1).Timestamp;
    zposArr = [zposArr round(datainfo.ScanFrames(i).ZPos*10)/10];

    % Later Additions (might not be present)
    datainfo.ScanFrames(i).DiodeReadout = nan;
    datainfo.ScanFrames(i).LaserEnergy = nan;
    try
        datainfo.ScanFrames(i).LaserEnergy = fr.getLASERENERGY.getDoubleValue;
    end
    try
        datainfo.ScanFrames(i).DiodeReadout = fr.getDIODEREADOUT.getDoubleValue;
    end
    
    ru = fr.getRUN;
    if ru == 0, ru = 1; end;            % if no RUNs, then it will be zero
    datainfo.ScanFrames(i).Run = ru;
    if (ru > maxRun) maxRun = ru; end
    datainfo.ScanFrames(i).Repetition = fr.getREPETITION;
    if (fr.getREPETITION > maxRep) maxRep = fr.getREPETITION; end
end;
datainfo.RepNum = maxRep;
datainfo.RunNum = maxRun; 
datainfo.ZPositions = unique(zposArr');
datainfo.AverageTemperature = tempavg / length(fa);
frmaxidoffset = fa(end).getIDOffset;

%% construct id and timestamp matrix
waitbar(0,wbar,'Verifying integrity...' ) ;

nzp = numel(datainfo.ZPositions);
numshots = datainfo.MeasurementDesc.NumFrames;
% Legacy support for datasets without Number of Frames
if (numshots == 0)
   numshots = round(numel(datainfo.ScanFrames)./(datainfo.RunNum*numel(datainfo.ZPositions)*datainfo.RepNum*numel(datainfo.Wavelengths)));
   datainfo.MeasurementDesc.NumFrames = numshots;
end
nestframes = numshots*numel(datainfo.Wavelengths)*datainfo.RepNum*nzp*datainfo.RunNum;
% compare actual number of frames with expected
if nestframes < numel(datainfo.ScanFrames),
   warning(['Corrupt data description - actual number of frames ' num2str(numel(datainfo.ScanFrames)) ' does not match expectations (' num2str(nestframes) ')']);
   nzpdiff = round((numel(datainfo.ScanFrames) - nestframes)/numel(datainfo.Wavelengths)/numshots);
   nzp = nzp + nzpdiff;
end
idm = nan(datainfo.RunNum,nzp,datainfo.RepNum,numel(datainfo.Wavelengths),numshots,'single');
tsm = nan(datainfo.RunNum,nzp,datainfo.RepNum,numel(datainfo.Wavelengths),numshots,'double');
lpm = nan(datainfo.RunNum,nzp,datainfo.RepNum,numel(datainfo.Wavelengths),numshots,'double');
idoffset = zeros(frmaxidoffset,1,'uint32');
sh = 1;
for i = 1:numel(datainfo.ScanFrames)
    waitbar(i/numel(datainfo.ScanFrames),wbar ) ;
    idoffset(datainfo.ScanFrames(i).IDOffset+1) = i;  % map id offset to frame number
    ru = datainfo.ScanFrames(i).Run;
    zp = find(datainfo.ZPositions == round(datainfo.ScanFrames(i).ZPos*10)/10);
    zp = zp(end);       % take the last of the found zpositions
    re = datainfo.ScanFrames(i).Repetition;
    wl = find(datainfo.Wavelengths == datainfo.ScanFrames(i).Wavelength);
    sh = 1;
    while (~isnan(idm(ru,zp,re,wl,sh))) 
        sh = sh + 1; 
        % catch if number of shots gets larger than anticipated
        % usually happens if stage unsuccessfully moves (which is a bug if
        % the tolarnce value is too low)
        % --> we create a duplicate the current zposition
        % MIGHT LEAD TO WEIRD BEHAVIOR IN OTHER CASES
        if sh > numshots
            if zp == numel(datainfo.ZPositions) || datainfo.ZPositions(zp) ~= datainfo.ZPositions(zp+1),
                datainfo.ZPositions = [datainfo.ZPositions(1:zp);datainfo.ZPositions(zp);datainfo.ZPositions(zp+1:end)];
                warning(['Corrupt data description in ' datainfo.Name ' (' num2str(ru) ',' num2str(zp) ',' num2str(re) ',' num2str(wl) ',' num2str(sh) ') - creating new virtual zpos (' num2str(datainfo.ZPositions(zp)) ')' ]);
            end
            sh = 1;
            zp = zp + 1;
            break
        end
    end
    if sh > numshots continue; end
    datainfo.ScanFrames(i).ZNum = zp;
    idm(ru,zp,re,wl,sh) = i;
    tsm(ru,zp,re,wl,sh) = datainfo.ScanFrames(i).RelTime;
    lpm(ru,zp,re,wl,sh) = datainfo.ScanFrames(i).LaserEnergy;
end
datainfo.ScanStructure = idm;
datainfo.RelTime = tsm;
datainfo.LaserEnergy = lpm;
datainfo.Timestamps = squeeze(datainfo.RelTime(:,1,1,1,1));
clear ru zp re wl sh u idm tsm;


%% recon nodes

waitbar(0,wbar,'Reading recon nodes...' ) ;
datainfo.ReconNode = struct([]);
rna = dm.getReconNodes.getDataModelNewReconstructionNodeArray;

% find total number of frames
rfrcount = 0;
rntotc = 0;
for rn = 1:length(rna)
    rfrcount = rfrcount + length(rna(rn).getReconFrames.getDataModelReconstructionFrameArray);
end
    
% read all Recons
for rn = 1:length(rna)
    rnode = rna(rn);
    datainfo.ReconNode(rn).Name = char(rnode.getName);
    datainfo.ReconNode(rn).Comment = char(rnode.getComment);
    datainfo.ReconNode(rn).GUID = char(rnode.getGUID);
    datainfo.ReconNode(rn).Method = char(rnode.getMethod);
    datainfo.ReconNode(rn).Resolution = rnode.getResolution;
%    datainfo.ReconNode(rn).ZResolution = rnode.getZResolution;
    datainfo.ReconNode(rn).ZResolution = rnode.getNumOfFrames;
    datainfo.ReconNode(rn).Projections = rnode.getProjections; 
    datainfo.ReconNode(rn).ROI = double(rnode.getRoi);
    datainfo.ReconNode(rn).ZPositions = [];
    datainfo.ReconNode(rn).ZNum = [];
    datainfo.ReconNode(rn).Runs = [];
    datainfo.ReconNode(rn).Repetitions = [];
    datainfo.ReconNode(rn).Wavelengths = [];
%     datainfo.ReconNode(rn).TrimSOS = rnode.getTrimSpeedOfSound;
    rfra = rnode.getReconFrames.getDataModelReconstructionFrameArray;
    % start ID matrix with same size as scan
%     ridm = nan(size(datainfo.ScanStructure),'single');
    clear ridm ridc
    sv = size(datainfo.ScanStructure);
    if numel(sv) < 4, sv = [sv ones(1,4-numel(sv))]; end;
    ridc = zeros(sv(1:4),'uint16');
    idlookup = zeros(length(rfra)*datainfo.ReconNode(rn).ZResolution,1,'uint16');
    for rfr = 1:length(rfra)
        rntotc = rntotc + 1;
        waitbar(rntotc/rfrcount,wbar ) ;
        datainfo.ReconNode(rn).Frames(rfr).IDOffset = rfra(rfr).getIDOffset;
        datainfo.ReconNode(rn).Frames(rfr).ScanFrames = str2num(rfra(rfr).getScanFramesID.getStringValue);     
        % use wavelength and zpos from first frame
        frid = idoffset(datainfo.ReconNode(rn).Frames(rfr).ScanFrames(1)+1);
        datainfo.ReconNode(rn).Frames(rfr).Wavelength = datainfo.ScanFrames(frid).Wavelength;
        datainfo.ReconNode(rn).Frames(rfr).ZPos = datainfo.ScanFrames(frid).ZPos;
        datainfo.ReconNode(rn).Frames(rfr).ZNum = datainfo.ScanFrames(frid).ZNum;
        datainfo.ReconNode(rn).Frames(rfr).Run = datainfo.ScanFrames(frid).Run;
        datainfo.ReconNode(rn).Frames(rfr).Repetition = datainfo.ScanFrames(frid).Repetition;
        datainfo.ReconNode(rn).Frames(rfr).Timestamp = datainfo.ScanFrames(frid).Timestamp;
        datainfo.ReconNode(rn).Frames(rfr).RelTime = datainfo.ScanFrames(frid).RelTime;
        
        % get all power values
        for j = 1:length(datainfo.ReconNode(rn).Frames(rfr).ScanFrames)
            datainfo.ReconNode(rn).Frames(rfr).LaserEnergy(j) = datainfo.ScanFrames(frid).LaserEnergy;
        end
        
        % assign id to correct entry in structure
        ru = find(datainfo.ReconNode(rn).Runs == datainfo.ReconNode(rn).Frames(rfr).Run);
        if isempty(ru), 
            ru = numel(datainfo.ReconNode(rn).Runs)+1; 
            datainfo.ReconNode(rn).Runs = [datainfo.ReconNode(rn).Runs, datainfo.ReconNode(rn).Frames(rfr).Run]; 
        end
        zp = find(datainfo.ReconNode(rn).ZNum == datainfo.ReconNode(rn).Frames(rfr).ZNum);
        if isempty(zp), 
            zp = numel(datainfo.ReconNode(rn).ZNum)+1; 
            datainfo.ReconNode(rn).ZNum = [datainfo.ReconNode(rn).ZNum datainfo.ReconNode(rn).Frames(rfr).ZNum]; 
        end
        re = find(datainfo.ReconNode(rn).Repetitions == datainfo.ReconNode(rn).Frames(rfr).Repetition);
        if isempty(re), 
            re = numel(datainfo.ReconNode(rn).Repetitions)+1; 
            datainfo.ReconNode(rn).Repetitions = [datainfo.ReconNode(rn).Repetitions datainfo.ReconNode(rn).Frames(rfr).Repetition]; 
        end
        wl = find(datainfo.ReconNode(rn).Wavelengths == datainfo.ReconNode(rn).Frames(rfr).Wavelength);
        if isempty(wl), 
            wl = numel(datainfo.ReconNode(rn).Wavelengths)+1; 
            datainfo.ReconNode(rn).Wavelengths = [datainfo.ReconNode(rn).Wavelengths datainfo.ReconNode(rn).Frames(rfr).Wavelength]; 
        end
        % increment shot counter (because we lack another intrinsic
        % numbering in the XML
        ridc(ru,zp,re,wl) = ridc(ru,zp,re,wl) + 1;
        sh = ridc(ru,zp,re,wl);
        if sh > numshots
            warning(['Corrupt data description (' num2str(ru) ',' num2str(zp) ',' num2str(re) ',' num2str(wl) ',' num2str(sh) ')']);
            continue;
        end
        ridm(ru,zp,re,wl,sh) = rfr;
        tsm(ru,zp,re,wl,sh) = datainfo.ReconNode(rn).Frames(rfr).RelTime;
        idlookup(datainfo.ReconNode(rn).Frames(rfr).IDOffset+1) = rfr;    % reverse lookup to find when loading MSP in 3D cases, where IDOffset != Frame Number
        clear ru zp re wl sh;
    end
    % save ReconStructure as matrix
    datainfo.ReconNode(rn).ZPositions = datainfo.ZPositions(datainfo.ReconNode(rn).ZNum);
    datainfo.ReconNode(rn).ReconStructure = ridm;
    datainfo.ReconNode(rn).RelTime = tsm;
    datainfo.ReconNode(rn).Timestamps = squeeze(datainfo.RelTime(:,1,1,1,1));
    datainfo.ReconNode(rn).IDLookup = idlookup;
    clear ridm;
    
end

%% MSP nodes
waitbar(0,wbar,'Reading MSP nodes...' ) ;
datainfo.MSPNode = struct([]);
mna = dm.getMspNodes.getDataModelNewMspNodeArray;

% Find total number of slices
mfrcount = 0;
mntotc = 0;
for mn = 1:length(mna)
    mfrcount = mfrcount + length(mna(mn).getSlicesFrame.getMspSliceFrameArray);
end




for mn = 1:length(mna)
    waitbar(mn/length(mna),wbar ) ;
    mnode = mna(mn);
    datainfo.MSPNode(mn).Name = char(mnode.getName);
    datainfo.MSPNode(mn).Comment = char(mnode.getComment);
    datainfo.MSPNode(mn).Method = char(mnode.getMethod);
    datainfo.MSPNode(mn).GUID = char(mnode.getGUID);
    datainfo.MSPNode(mn).ReconGUID = char(mnode.getReconGUID);
    % Related Wavelengths
    wla = mnode.getRelatedWavelengths.getWavelengthArray;
    for wl = 1:numel(wla)
        datainfo.MSPNode(mn).Wavelengths(wl) = wla(wl).getIntValue;
    end
    % Input Spectra
    datainfo.MSPNode(mn).InputSpectra = cell(mnode.getInputSpectra.getStringArray);
    % Recon Info
    rn = 1;
    while ~strcmp(datainfo.ReconNode(rn).GUID,datainfo.MSPNode(mn).ReconGUID)
        rn = rn + 1;
        if (rn > numel(datainfo.ReconNode)) rn = 0; break; end;
    end
    datainfo.MSPNode(mn).ReconNodeID = rn;
    if (rn > 0)
        datainfo.MSPNode(mn).ReconMethod = datainfo.ReconNode(rn).Method;
        datainfo.MSPNode(mn).Resolution = datainfo.ReconNode(rn).Resolution;
        datainfo.MSPNode(mn).Projections = datainfo.ReconNode(rn).Projections;
        datainfo.MSPNode(mn).ROI = datainfo.ReconNode(rn).ROI;
    else
        warning(['Recon ' datainfo.MSPNode(mn).ReconGUID ' not found.\n']);
    end
        
    % SliceFrames
    datainfo.MSPNode(mn).ZPositions = [];  % track z positions
    datainfo.MSPNode(mn).ZNum = [];  % track z positions
    mfra = mnode.getSlicesFrame.getMspSliceFrameArray;
    clear mstr slref;
    for mfr = 1:numel(mfra)
        mntotc = mntotc + 1;
        waitbar(mntotc/mfrcount,wbar ) ;
        
        datainfo.MSPNode(mn).Slices(mfr).HasErrors = logical(mfra(mfr).getHasErrors);
        % Recon Info
        if rn > 0,
            datainfo.MSPNode(mn).Slices(mfr).ReconFrames = str2num(mfra(mfr).getReconFrameID.getStringValue);
            % reverse lookup to find recon frame in 3D cases, where IDOffset != Frame Number
            datainfo.MSPNode(mn).Slices(mfr).ReconFrameNumbers = datainfo.ReconNode(rn).IDLookup(datainfo.MSPNode(mn).Slices(mfr).ReconFrames(:)+1);
            rframeid = datainfo.ReconNode(rn).IDLookup(datainfo.MSPNode(mn).Slices(mfr).ReconFrames(1)+1);
            datainfo.MSPNode(mn).Slices(mfr).ZPos = datainfo.ScanFrames(rframeid).ZPos;
            datainfo.MSPNode(mn).Slices(mfr).ZNum = datainfo.ScanFrames(rframeid).ZNum;
            datainfo.MSPNode(mn).Slices(mfr).Run = datainfo.ScanFrames(rframeid).Run;
            datainfo.MSPNode(mn).Slices(mfr).Repetition = datainfo.ScanFrames(rframeid).Repetition;
            datainfo.MSPNode(mn).Slices(mfr).RelTime = datainfo.ScanFrames(rframeid).RelTime;
            zp = find(datainfo.MSPNode(mn).ZNum ==  datainfo.MSPNode(mn).Slices(mfr).ZNum);
            % if z-position was not yet added, add it now and keep the index of
            % the added z-positions for the structure indexing
            if (isempty(zp)) 
                datainfo.MSPNode(mn).ZNum = [datainfo.MSPNode(mn).ZNum datainfo.MSPNode(mn).Slices(mfr).ZNum];
                zp = numel(datainfo.MSPNode(mn).ZNum);
            end
        end
        % Component Info
        ca = mfra(mfr).getComponentList.getMspSliceComponentArray;
        for c = 1:numel(ca)
           datainfo.MSPNode(mn).Slices(mfr).Components(c).IDOffset = ca(c).getIDOffset;
           datainfo.MSPNode(mn).Slices(mfr).Components(c).ComponName = char(ca(c).getComponName);
           datainfo.MSPNode(mn).Slices(mfr).Components(c).ComponRef = char(ca(c).getComponRef);
           datainfo.MSPNode(mn).Slices(mfr).Components(c).Invert = logical(ca(c).getInvert);
           datainfo.MSPNode(mn).Slices(mfr).Components(c).SpectrumName = char(ca(c).getSpectraName);
           mstr(datainfo.MSPNode(mn).Slices(mfr).Run,zp,datainfo.MSPNode(mn).Slices(mfr).Repetition,c) = ca(c).getIDOffset+1;
           slref(ca(c).getIDOffset+1) = mfr;
        end
    end
    datainfo.MSPNode(mn).ZPositions = datainfo.ZPositions(datainfo.MSPNode(mn).ZNum);
    datainfo.MSPNode(mn).Structure = mstr;
    datainfo.MSPNode(mn).SliceIndex = slref;
end

%%
close(wbar);


return


end











