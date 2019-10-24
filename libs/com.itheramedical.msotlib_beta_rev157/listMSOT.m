function datainfoarr  = listMSOT( varargin )
% LISTMSOT  List .msot META information from Study Folder
%   datainfo = LISTMSOT()      List scans in current directory
%   datainfo = LISTMSOT(dir)   Fist scans in the specified directory
% 
% This function lists the META information on the selected study folder
%
% Required Parameters:
%   <none>
%
% Optional:
%   1: dirname:    path of the study to load.
%
% Return Values:
%   1: datainfoarr: a structure containing the meta information, [] on error 
%
% Example:
%   datainfo = loadMSOT('Scan_1\Scan_1.msot');


%% parameters
if ( nargin == 0 )
   studydir = '.'; 
elseif( nargin >= 1 )
   studydir = varargin{1} ;
else
   datainfo = [];
   error( 'Not enough input arguments' ) ;
end ;



%% file list
wbar = waitbar(0,'Initialising...');
cfiles = rdir([studydir '\**\' '*.mso*']);

for ind = 1:numel(cfiles)
    xml_filename = cfiles(ind).name;
    bsl = strfind(xml_filename,'\');
    dirname = xml_filename(1:bsl(end)-1);
    
    waitbar(ind/numel(cfiles),wbar,'Parsing XML...');
    try
        patdoc = javaMethod( 'parse', 'com.itheramedical.msotbeans.DataModelMsotProjectDocument$Factory', java.io.File( xml_filename ) ) ;
    catch ex
        close(wbar)
        error( ['Selected XML file could not be loaded. Check well-formedness.\n' ex.message]) ;
    end ;
    dm = patdoc.getDataModelMsotProject() ;

    datainfoarr(ind).DirName = char(dm.getFriendlyName);
    datainfoarr(ind).FileName = xml_filename;
    datainfoarr(ind).CreationTimeTxt = char(dm.getCreationTime.toString);
    datainfoarr(ind).CreationTime = datenum(datainfoarr(ind).CreationTimeTxt,'yyyy-mm-ddTHH:MM:SS.FFF');
    datainfoarr(ind).Name = char(dm.getScanNode.getName);
    datainfoarr(ind).RealPath = dirname;

    hw = dm.getHARDWAREDESC;
    datainfoarr(ind).HWDesc.TransducerType = char(hw.getTRANSDUCER);
    datainfoarr(ind).is3D = strcmp(datainfoarr(ind).HWDesc.TransducerType,'msot3');
    datainfoarr(ind).is2D = strcmp(datainfoarr(ind).HWDesc.TransducerType,'msot2');
    if (datainfoarr(ind).is3D),
        datainfoarr(ind).TType = '3D';
    else
        datainfoarr(ind).TType = '2D';
    end

    
    rna = dm.getReconNodes.getDataModelNewReconstructionNodeArray;
    datainfoarr(ind).ReconNode = [];
    for rn = 1:length(rna)
       rnode = rna(rn);
       datainfoarr(ind).ReconNode(rn).Name = char(rnode.getName);
       datainfoarr(ind).ReconNode(rn).Comment = char(rnode.getComment);
       datainfoarr(ind).ReconNode(rn).GUID = char(rnode.getGUID);
       datainfoarr(ind).ReconNode(rn).Method = char(rnode.getMethod);
       datainfoarr(ind).ReconNode(rn).Resolution = rnode.getResolution;
       datainfoarr(ind).ReconNode(rn).Projections = rnode.getProjections;
       datainfoarr(ind).ReconNode(rn).ROI = double(rnode.getRoi);
    end
    
    
    mna = dm.getMspNodes.getDataModelNewMspNodeArray;
    for mn = 1:length(mna)
        mnode = mna(mn);
        datainfoarr(ind).MSPNode(mn).Name = char(mnode.getName);
        datainfoarr(ind).MSPNode(mn).Comment = char(mnode.getComment);
        datainfoarr(ind).MSPNode(mn).Method = char(mnode.getMethod);
        datainfoarr(ind).MSPNode(mn).GUID = char(mnode.getGUID);
        datainfoarr(ind).MSPNode(mn).ReconGUID = char(mnode.getReconGUID);
        % Related Wavelengths
        wla = mnode.getRelatedWavelengths.getWavelengthArray;
        for wl = 1:numel(wla)
            datainfoarr(ind).MSPNode(mn).Wavelengths(wl) = wla(wl).getIntValue;
        end
        % Input Spectra
        datainfoarr(ind).MSPNode(mn).InputSpectra = cell(mnode.getInputSpectra.getStringArray);
        % Recon Info
        rn = 1;
        while ~strcmp(datainfoarr(ind).ReconNode(rn).GUID,datainfoarr(ind).MSPNode(mn).ReconGUID)
            rn = rn + 1;
            if (rn > numel(datainfoarr(ind).ReconNode)) rn = 0; break; end;
        end
        datainfoarr(ind).MSPNode(mn).ReconNodeID = rn;
        if (rn > 0)
            datainfoarr(ind).MSPNode(mn).ReconMethod = datainfoarr(ind).ReconNode(rn).Method;
            datainfoarr(ind).MSPNode(mn).Resolution = datainfoarr(ind).ReconNode(rn).Resolution;
            datainfoarr(ind).MSPNode(mn).Projections = datainfoarr(ind).ReconNode(rn).Projections;
            datainfoarr(ind).MSPNode(mn).ROI = datainfoarr(ind).ReconNode(rn).ROI;
        else
            warning(['Recon ' datainfoarr(ind).MSPNode(mn).ReconGUID ' not found.\n']);
        end
    end
     
    
end
close(wbar);


%%
[tmp ind] = sort(arrayfun(@(x) x.CreationTime, datainfoarr));
datainfoarr = datainfoarr(ind);

for ind = 1:numel(datainfoarr)
    fprintf('%s\t\t%s\t%s\t%s\n',datainfoarr(ind).DirName,datainfoarr(ind).CreationTimeTxt,datainfoarr(ind).TType,datainfoarr(ind).Name);
end

return


