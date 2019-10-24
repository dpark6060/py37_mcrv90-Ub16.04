function [ sigMat, datainfo ] = loadData( varargin )
%% Data Loader - loads acoustic signal data from a binary file with xml header.
%
% @param filename String, name of the file to load.
% @param format String, data format: 'mat' or 'xml'.
% @param progressBar, javax.swing.JProgressBar to set the loading progress
%                     in the GUI Dataloader dialog
% @param partialLoadingFlag, disable partial data loading. Enabled by default otherwise
%
% @return sigMat, matrix containing all projections. Size of sigMat is 
%                 (m, n, z), where m is the number of timesamples per
%                 projection, n the number of projections and z the number
%                 of measurements in the file (i.e. slices, vertical steps).
% @return datainfo, a structure containing the measurement information as
%                   stored in the datafile.

global root ;
    
if ( nargin == 0 )
    
    path = getappdata( 0, 'last_path' ) ;
    if( isempty( path ) )
        path = cd ;
    end ;
    fc = javax.swing.JFileChooser( java.io.File( path ) ) ;
    fc.setAcceptAllFileFilterUsed(0) ;
    fc.addChoosableFileFilter( javax.swing.filechooser.FileNameExtensionFilter(  'XML Header',  'xml' ) ) ;
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
    progressBar = [] ;
    partialLoadingFlag = 1 ;
    
elseif( nargin == 2 )
    filename = varargin{1} ;
    format = varargin{2} ;
    progressBar = [] ;
    partialLoadingFlag = 1 ;
elseif( nargin == 3 )
    filename = varargin{1} ;
    format = varargin{2} ;
    progressBar = varargin{3} ;
    partialLoadingFlag = 1 ;
elseif( nargin == 4 )
    filename = varargin{1} ;
    format = varargin{2} ;
    progressBar = varargin{3} ;
    partialLoadingFlag = varargin{4} ;
else
    error( 'Not enough input arguments' ) ;
end ;

if( strmatch( format, 'mat', 'exact' ) )
    
    vars = whos('-file', filename) ;
    found = 0 ;
    found_info = 0 ;
    for ii = 1:length( vars )
        if( strmatch( vars(ii).name, 'sigMat', 'exact' ) )
            found = 1 ;
        end ;
        if( strmatch( vars(ii).name, 'datainfo', 'exact' ) )
            found_info = 1 ;
        end ;
    end ;
    if( ~found )
        javax.swing.JOptionPane.showMessageDialog( [], [ 'No signal matrix found in file ' filename '.' ] ) ;
    end ;
    if( ~found_info )
%         javax.swing.JOptionPane.showMessageDialog( [], [ 'No datainfo found in file ' filename '.' ] ) ;
        datainfo = struct( 'isXML', 0 ) ;
    else
        data = load( filename, 'datainfo' ) ;
        datainfo = data.datainfo ;
    end ;
    
    data = load( filename, 'sigMat' ) ;
    sigMat = data.sigMat ;
    return ;
end ;

dlg = [] ;
if( isempty( progressBar ) )
    if( isempty( getappdata( 0, 'loading_dialog' ) ) )
        dlg = com.helmholtz.pat.PATDaqProgressDialog( [], 0 ) ;
        dlg.getContentPane().getComponent(1).setText( 'Hide' ) ;
        dlg.setTitle( 'Load data progress' ) ;
        setappdata( 0, 'loading_dialog', dlg ) ;
    else
        dlg = getappdata( 0, 'loading_dialog' ) ;
    end ;
    progressBar = dlg.getContentPane().getComponent(0) ;
    dlg.setVisible(1) ;
end ;
progressBar.setValue( 0 ) ;

fstr = java.lang.String( filename ) ;
xml_filename = '' ;
if( fstr.endsWith( '.xml' ) )
    xml_filename = filename ;
    filename = char( fstr.substring( 0, fstr.lastIndexOf( '.xml' ) ).concat( '.bin' ) ) ;
end ;

if( ~java.io.File( filename ).exists() )
    error( [ 'File ' char( filename ) ' does not exist!' ] ) ;
end ;

%% XML header format
FID = fopen(filename, 'r', 'ieee-le') ;
try
    tic
    patdoc = javaMethod( 'parse', 'com.helmholtz.pat.patbeans.PATROOTDocument$Factory', java.io.File( xml_filename ) ) ;
    toc
catch ex
    error( 'Selected XML file could not be loaded. Check well-formedness.' ) ;
end ;
if( ~patdoc.validate() )
    error( 'XML document not valid. Validate against Header_Schema.xsd' ) ;
end ;
root = patdoc.getPATROOT() ;
%% declaration
datainfo = struct(  'admin_data'        ,  struct(),   ...
                    'hardware_desc'     ,  struct(),   ...
                    'measurement_desc'  ,  struct(),   ...
                    'comment'           ,  char( root.getCOMMENT() ), ...
                    'isXML'             ,  1 ) ;

progressBar.setIndeterminate(1) ;
progressBar.setString( 'Parsing XML, please wait.' ) ;
                
%% admin data
tmp = root.getADMINDATA().getDATE() ;
datainfo.admin_data.date = [ num2str( tmp.getDAY(), '%02d' ) '.' num2str( tmp.getMONTH(), '%02d' ) '.' num2str( tmp.getYEAR(), '%02d' ) ] ;
datainfo.admin_data.target = char( root.getADMINDATA().getMEASUREMENTTARGET() ) ;
tma = root.getADMINDATA().getTEAMMEMBERArray() ;
owner_str = '' ;
for ii = 1 : length( tma )
    owner_str = [ owner_str char( tma(ii) ) ] ; %#ok<AGROW>
    if( ii ~= length( tma ) )
        owner_str = [ owner_str ', ' ] ; %#ok<AGROW>
    end ;
end ;
datainfo.admin_data.owner = owner_str ;

%% hardware description
datainfo.hardware_desc = readHWDesc( root.getHARDWAREDESC() ) ;

%% measurement description
datainfo.measurement_desc = readMeasurementDesc( root.getMEASUREMENTDESC(), datainfo ) ;

progressBar.setIndeterminate(0) ;
progressBar.setString( [] ) ;
 
%% read binary file
max_wavelengths = 1 ;
if( isfield( datainfo.measurement_desc, 'wavelengths' ) )
    max_wavelengths = size( datainfo.measurement_desc.wavelengths, 1 ) ;
end ;
samples = datainfo.measurement_desc.rec_length ;
proj = datainfo.measurement_desc.projections.num_proj ; % number of projections per frame
if( isfield( datainfo.measurement_desc.projections, 'equal' ) )
    slices = length( datainfo.measurement_desc.projections.equal ) ;
    frames = datainfo.measurement_desc.projections.equal ;
elseif( isfield( datainfo.measurement_desc.projections, 'msot_frame' ) )
    slices = length( datainfo.measurement_desc.projections.msot_frame ) ;
    frames = datainfo.measurement_desc.projections.msot_frame ;
    max_wavelengths = 1 ;
else
    % the PROJECTION choice is treated as a PROJECTION frame containing
    % one projection
    slices = length( datainfo.measurement_desc.projections.projection ) ;
    frames = datainfo.measurement_desc.projections.projection ;
end ;

% if the data to be loaded exceeds 1 GB of RAM, do partial data loading
% Implemented only for MSOT-FRAME definition
if( partialLoadingFlag && (samples*proj*slices*max_wavelengths *8/1e9 > .5 || ...
                    ~javax.swing.JOptionPane.showConfirmDialog( [], 'Partial data loading?', 'Data loader', javax.swing.JOptionPane.YES_NO_OPTION ) ) )
    
    num_values = samples*proj ;
    
    pldDlg = com.helmholtz.pat.PATLoadDataDlg.getInstance() ;
    for ii = 1 : slices
        pldDlg.addFrame( frames{ii}.zpos, frames{ii}.actual_wavelength*1e9, ii ) ;
    end ;
    
    selection_valid = 0 ;
    while( ~selection_valid )
        pldDlg.setVisible(1) ;
        % user input @modal dialog
        selection = double( pldDlg.getSelection() ) ;
        if( length(selection) * num_values *8/1e9 > .5 )
            yesno = javax.swing.JOptionPane.showConfirmDialog( [], ['Selected data will occupy ' num2str( length(selection) * num_values *8/1e9 ) ' GB in memory. Proceed?'], 'Confirm selection', javax.swing.JOptionPane.YES_NO_OPTION ) ;
            if( ~yesno )
                % proceed
                selection_valid = 1 ;
            end ;
        else
            selection_valid = 1 ;
        end ;
    end ;
    
    progressBar.setMaximum( length( selection ) ) ;
        
    sigMat = zeros( samples, proj, length( selection ) ) ;
    
    if( pldDlg.cancelFlag )
        if( ~isempty( getappdata( 0, 'loading_dialog' ) ) ) 
            dlg.setVisible(0) ;
        end ;
        return ;
    end ;
    
    frame_counter = 1 ;
    for ii = selection'
        
        if( ftell( FID ) ~= 2* num_values * (ii-1) )
            fseek( FID, 2*num_values*(ii-1), 'bof' ) ;
        end ;
        sigMat( :, :, frame_counter ) = reshape( fread( FID, num_values, 'uint16' ), samples, proj ) ;
        frame_counter = frame_counter + 1 ;
        progressBar.setValue( progressBar.getValue() + 1 ) ;
        
    end ;
    datainfo.measurement_desc.projections.msot_frame = datainfo.measurement_desc.projections.msot_frame(selection) ;
else
    
    progressBar.setMaximum( length( frames ) * max_wavelengths ) ;
    
    sigMat = zeros( samples, proj, slices, max_wavelengths ) ;
    
    if( ~isempty( strmatch( datainfo.measurement_desc.sequence, 'by-wavelength', 'exact' ) ) )
        % datastream order:
        % same frame in multiple (sequential) wavelengths
        for ii = 1 : length( frames )
            curr_frame = frames{ ii, 1 } ;

            if( isfield( curr_frame, 'wavelengths' ) )
                curr_wavelengths = curr_frame.wavelengths ;
            else
                curr_wavelengths = 1 ;
            end ;

            num_proj = proj ;
            if( isfield( curr_frame, 'num_proj' ) )
                num_proj = curr_frame.num_proj ;
            end ;
            
            for jj = 1 : size( curr_wavelengths, 1 )

%               SMorscher, 06/08/2012: * Handle corrupt binary files
                tmp = fread( FID, samples*num_proj, 'uint16' );
                try
                    % this may fail if crash occurred during acquisition
                    sigMat( :, :, ii, curr_wavelengths( jj, 1 ) ) = reshape( tmp, samples, num_proj ) ;
                catch err
                    % print warning, but leave the data as initialised
                    fprintf( 'WARNING: Frame data incomplete, putting to zero. (%i bytes of %i total)\n   MATLAB message: %s\n', 2*size( tmp, 1 ), samples * num_proj, err.message ) ;
                    break;
                end;
                
                progressBar.setValue( progressBar.getValue() + 1 ) ;

            end ;
            clear tmp

        end ;

    else
        % datastream order:
        % same wavelength for multiple (sequential) frames
        for ii = 1 : max_wavelengths

            for jj = 1 : length( frames )
                curr_frame = frames{ jj, 1 } ;

                if( isfield( curr_frame, 'wavelengths' ) )
                    wlen_vec = curr_frame.wavelengths( :, 1 ) ;
                else
                    wlen_vec = 1 ;
                end ;

                if( ii > length( wlen_vec ) )
                    continue ;
                end ;

                num_proj = proj ;
                if( isfield( curr_frame, 'num_proj' ) )
                    num_proj = curr_frame.num_proj ;
                end ;

                sigMat( :, :, jj, wlen_vec( ii ) ) = reshape( fread( FID, samples*num_proj, 'uint16' ), samples, num_proj ) ;

                progressBar.setValue( progressBar.getValue() + 1 ) ;

            end ;
        end ;

    end ;
end ;

for ii = 1 : size( sigMat, 3 )
    sigMat( :, :, ii ) = sigMat( :, :, ii ) / datainfo.measurement_desc.projections.msot_frame{ii}.actual_power ;
end ;

if( ~isempty( dlg ) ) 
    dlg.setVisible(0) ;
end ;

clear global root ;
clear functions ;

fclose( FID ) ;

end

function ms = readMeasurementDesc( msbean, datainfo )
    % required
    ms = struct( 'rec_length'   ,   0, ...
                 'sequence'     ,  '', ...
                 'coord_def'    ,  '', ...
                 'projections'  ,  struct() ) ;
    
    ms.rec_length  = msbean.getRECORDEDLENGTH() ;
    ms.sequence    = char( msbean.getSEQUENCE().toString() ) ;
    
    linesep = char( java.lang.System.getProperty( 'line.separator' ) ) ;
    unit_mod = msbean.getCOORDINATEDEFINITION().getAXISArray(0).getUnitModifier() ;
    if( unit_mod.equals( 'none' ) )
        unit_mod = '' ;
    end ;
    str = [ char( msbean.getCOORDINATEDEFINITION().getAXISArray(0).getStringValue() ) ...
        ' in ' ...
        char( unit_mod ) ...
        char( msbean.getCOORDINATEDEFINITION().getAXISArray(0).getUnit() ) ...
        linesep ] ;
    
    unit_mod = msbean.getCOORDINATEDEFINITION().getAXISArray(1).getUnitModifier() ;
    if( unit_mod.equals( 'none' ) )
        unit_mod = '' ;
    end ;
    str = [ str ...
        char( msbean.getCOORDINATEDEFINITION().getAXISArray(1).getStringValue() ) ...
        ' in ' ...
        char( unit_mod ) ...
        char( msbean.getCOORDINATEDEFINITION().getAXISArray(1).getUnit() ) ...
        linesep ] ;
    
    unit_mod = msbean.getCOORDINATEDEFINITION().getAXISArray(2).getUnitModifier() ;
    if( unit_mod.equals( 'none' ) )
        unit_mod = '' ;
    end ;
    str = [ str ...
        char( msbean.getCOORDINATEDEFINITION().getAXISArray(2).getStringValue() ) ...
        ' in ' ...
        char( unit_mod ) ...
        char( msbean.getCOORDINATEDEFINITION().getAXISArray(2).getUnit() ) ...
        linesep ] ;
    
    ms.coord_def   = str ;
    
    ms.projections = getDetectorPositions( msbean.getPROJECTIONS(), datainfo ) ;
    
    % optional
    if( msbean.isSetTEMPERATURE() )
        ms.temperature = msbean.getTEMPERATURE().getDoubleValue() * getFactorFromENG( msbean.getTEMPERATURE().getUnitModifier() ) ;
    end ;
    if( msbean.isSetAVERAGESPERPROJECTION() )
        ms.averages = msbean.getAVERAGESPERPROJECTION() ;
    end ;
    if( msbean.isSetREPETITIONRATE() )
        ms.rep_rate = msbean.getREPETITIONRATE().getDoubleValue() * getFactorFromENG( msbean.getREPETITIONRATE().getUnitModifier() ) ; 
    end ;
    if( msbean.isSetVERTICALDEPTH() )
        ms.vert_depth = msbean.getVERTICALDEPTH().getDoubleValue() * getFactorFromENG( msbean.getVERTICALDEPTH().getUnitModifier() ) ; 
    end ;
    if( msbean.isSetINITIALENERGY() )
        ms.init_energy = msbean.getINITIALENERGY().getDoubleValue() * getFactorFromENG( msbean.getINITIALENERGY().getUnitModifier() ) ; 
    end ;
    if( msbean.isSetVERTICALSTEPS() )
        ms.vert_steps = msbean.getVERTICALSTEPS() ; 
    end ;
    if( msbean.isSetWAVELENGTHS() )
        ms.wavelengths = resolveWavelengthRef( java.lang.String( 'all' ) ) ;
    end ;
    
end 

function hwd = readHWDesc( hwbean )
    % required
    hwd = struct(  'fs',             0, ...
                   'geometry',      '', ...
                   'AD_range',       0, ...
                   'transducer',    '', ...
                   'power',          0, ...
                   'setup_type',    '' ) ;
    
    hwd.fs          = hwbean.getSAMPLINGFREQUENCY().getDoubleValue() * getFactorFromENG( hwbean.getSAMPLINGFREQUENCY().getUnitModifier() ) ;
    hwd.AD_range    = hwbean.getADRANGE().getDoubleValue() * getFactorFromENG( hwbean.getADRANGE().getUnitModifier() ) ;
    hwd.power       = hwbean.getPOWER().getDoubleValue() * getFactorFromENG( hwbean.getPOWER().getUnitModifier() ) ;
    
    hwd.geometry = char( hwbean.getGEOMETRY() ) ;
    hwd.transducer = char( hwbean.getTRANSDUCER() ) ;
    hwd.setup_type = char( hwbean.getSETUPTYPE() ) ;
    
    % optional
    if( hwbean.isSetAMPLIFICATION() )
        hwd.amplification = hwbean.getAMPLIFICATION() ;
    end ;
    if( hwbean.isSetVERTICALRANGEPHOTODIODE() )
        hwd.vert_range_diode = hwbean.getVERTICALRANGEPHOTODIODE().getDoubleValue() * getFactorFromENG( hwbean.getVERTICALRANGEPHOTODIODE().getUnitModifier() ) ;
    end ;
    if( hwbean.isSetFRAMEDESC() )
        hwd.framedesc = getDetectorPositions( hwbean.getFRAMEDESC(), [] ) ;        
    end ;
end 

function posdef = getDetectorPositions( parent, datainfo )
    
    posdef = struct() ;
    posdef.num_proj = -1 ;

    if( ( parent.instanceType.getShortJavaName.equals( 'FRAMEDESC' ) && parent.isSetEQUAL() ) || ( ~parent.instanceType.getShortJavaName.equals( 'FRAMEDESC' ) && parent.sizeOfEQUALArray() ) )
        % EQUAL choice
        equalSize = 1 ;
        if( ~parent.instanceType.getShortJavaName.equals( 'FRAMEDESC' ) )
            equalSize = parent.sizeOfEQUALArray() ;
        end ;
        posdef.equal = cell( equalSize, 1 ) ;
        for ii = 1 : equalSize
            if( ~parent.instanceType.getShortJavaName.equals( 'FRAMEDESC' ) )
                e = parent.getEQUALArray( ii-1 ) ;
            else
                e = parent.getEQUAL() ;
            end ;
            
            tmp = resolveAxisValue( e.getCONSTANTArray(0).getAxisRef(), e.getCONSTANTArray(0).getDoubleValue() ) ;
            eval( [ 'curr_equal.const_' tmp.name ' = ' num2str( tmp.value ) ' ; ' ] ) ;
            tmp = resolveAxisValue( e.getCONSTANTArray(1).getAxisRef(), e.getCONSTANTArray(1).getDoubleValue() ) ;
            eval( [ 'curr_equal.const_' tmp.name ' = ' num2str( tmp.value ) ' ; ' ] ) ;
            
            tmp = resolveAxisValue( e.getAxisRef(), 1 ) ;
            eval( [ 'curr_equal.var_' tmp.name ' = ' num2str( tmp.value ) ' * [ ' num2str( e.getSTART(), '%1.24f' ) ':' num2str( e.getSTEP(), '%1.24f' ) ':' num2str( e.getEND(), '%1.24f' ) ' ] ; ' ] ) ;
            
            if( posdef.num_proj < length( eval( [ 'curr_equal.var_' tmp.name ] ) ) )
                posdef.num_proj = length( eval( [ 'curr_equal.var_' tmp.name ] ) ) ;
            end ;
            
            curr_equal.wavelengths = resolveWavelengthRef( e.getWavelengthRef() ) ;
            curr_equal.num_proj = length( eval( [ 'curr_equal.var_' tmp.name ] ) ) ;
            curr_equal.number = e.getNumber() ;
            
            posdef.equal{ii, 1} = curr_equal ;
            
        end ;
        
    elseif( parent.sizeOfPROJECTIONArray() )
        % PROJECTION choice
        posdef.projection = cell( parent.sizeOfPROJECTIONArray(), 1 ) ;
        for ii = 1 : parent.sizeOfPROJECTIONArray()
            
            p = parent.getPROJECTIONArray( ii-1 ) ;
            
            tmp = resolveAxisValue( p.getVALUEArray(0).getAxisRef(), p.getVALUEArray(0).getDoubleValue() ) ;
            eval( [ 'curr_proj.' tmp.name ' = ' num2str( tmp.value ) ' ;' ] ) ;
            tmp = resolveAxisValue( p.getVALUEArray(1).getAxisRef(), p.getVALUEArray(1).getDoubleValue() ) ;
            eval( [ 'curr_proj.' tmp.name ' = ' num2str( tmp.value ) ' ;' ] ) ;
            tmp = resolveAxisValue( p.getVALUEArray(2).getAxisRef(), p.getVALUEArray(2).getDoubleValue() ) ;
            eval( [ 'curr_proj.' tmp.name ' = ' num2str( tmp.value ) ' ;' ] ) ;
            
            curr_proj.wavelengths = resolveWavelengthRef( p.getWavelengthRef() ) ;
            curr_proj.number = p.getNumber() ;
            
            posdef.projection{ii, 1} = curr_proj ;
            
        end ;
        posdef.num_proj = length( posdef.projection ) ;
        
    elseif( parent.sizeOfMSOTFRAMEArray() )
        % MSOT-FRAME choice        
        if( ~isfield( datainfo.hardware_desc, 'framedesc' ) )
            error( 'Missing FRAME-DESC for MSOT-FRAME.' ) ;
        end ;
        
        posdef.msot_frame = cell( parent.sizeOfMSOTFRAMEArray(), 1 ) ;
        for ii = 1 : parent.sizeOfMSOTFRAMEArray()
            
            mf = parent.getMSOTFRAMEArray( ii-1 ) ;
            
            tmp = resolveAxisValue( mf.getZPOS().getAxisRef(), mf.getZPOS().getDoubleValue() ) ;
            curr_frame.zpos = tmp.value ;
            curr_frame.actual_temp = mf.getTEMPERATURE().getDoubleValue() * getFactorFromENG( mf.getTEMPERATURE().getUnitModifier() ) ;
            curr_frame.actual_power = mf.getPOWER().getDoubleValue() * getFactorFromENG( mf.getPOWER().getUnitModifier() ) ;
            curr_frame.actual_wavelength = mf.getWAVELENGTH().getDoubleValue() * getFactorFromENG( mf.getWAVELENGTH().getUnitModifier() ) ;
            curr_frame.timestamp = mf.getTimestamp() * 1e-4 ; % - 2010 * 365 * 24 * 3600 * 1e3 ; % starting in the year 2010
            
            posdef.msot_frame{ii, 1} = curr_frame ;
            
            if( posdef.num_proj < 0 )
                posdef.num_proj = datainfo.hardware_desc.framedesc.num_proj ;
            end ;
            
        end ;
        
    end ;

end

function wavelengths = resolveWavelengthRef( refstr )
    global root ;
    persistent wlen_hmap ;
    if( isempty( wlen_hmap ) )
        hmSize = root.getMEASUREMENTDESC().getWAVELENGTHS().sizeOfWAVELENGTHArray() ;
        wlen_hmap = java.util.HashMap( hmSize ) ;
        wlens = root.getMEASUREMENTDESC().getWAVELENGTHS() ;
        for ii = 1 : hmSize
            wlen_hmap.put( wlens.getWAVELENGTHArray( ii-1 ).getNumber(), ...
                java.lang.Double( wlens.getWAVELENGTHArray( ii-1 ).getDoubleValue() * getFactorFromENG( wlens.getWAVELENGTHArray( ii-1 ).getUnitModifier() ) ) ) ;
        end ;
    end ;

    wavelengths = zeros( 1, 2 ) ;
    wlens = root.getMEASUREMENTDESC().getWAVELENGTHS() ;
    if( refstr.equals( 'all' ) )
        for ii = 1 : wlens.sizeOfWAVELENGTHArray()
            tmp = wlens.getWAVELENGTHArray( ii-1 ) ;
            wavelengths( ii, : ) = [ tmp.getNumber() tmp.getDoubleValue() * getFactorFromENG( tmp.getUnitModifier() ) ] ;
        end ;
        return ;
    end ;
    
    emptysp = java.lang.String(' ') ;
    curr_ind = refstr.indexOf( emptysp, 0 ) ;
    last_ind = 0 ;
    counter = 1 ;
    while( curr_ind > 0 )
        curr_substring = refstr.substring( last_ind, curr_ind ) ;
        if( ~wlen_hmap.containsKey( str2double( curr_substring ) ) )
            error( 'Referenced wavelength not defined! Check WAVELENGTHS.' ) ;
        end ;
        wavelengths( counter, : ) = [ str2double( curr_substring ) wlen_hmap.get( str2double( curr_substring ) ) ] ;
        counter = counter + 1 ;
        last_ind = curr_ind + 1 ;
        curr_ind = refstr.indexOf( emptysp, last_ind ) ;
    end ;
    curr_substring = refstr.substring( last_ind ) ;
    wavelengths( counter, : ) = [ str2double( curr_substring ) wlen_hmap.get( str2double( curr_substring ) ) ] ;
    
end

function axis_val = resolveAxisValue( ax_num, value )
    global root ;
    axs = root.getMEASUREMENTDESC().getCOORDINATEDEFINITION().getAXISArray() ;
    
    ax = [] ;
    for ii = 1 : length( axs )
        if( axs(ii).getNumber() == ax_num )
            ax = axs(ii) ;
            break ;
        end ;
    end ;
    if( isempty( ax ) )
        error( 'Invalid Axis-Ref! Check XML header.' ) ;
    end ;
    
    axis_val.name = char( ax.getStringValue() ) ;
    axis_val.value = value * getFactorFromENG( ax.getUnitModifier() ) ;
end

function factor = getFactorFromENG( str )
    if( str.equals( java.lang.String('Y') ) )
        factor = 1e24 ;
    elseif( str.equals( java.lang.String('Z') ) )
        factor = 1e21 ;
    elseif( str.equals( java.lang.String('E') ) )
        factor = 1e18 ;
    elseif( str.equals( java.lang.String('P') ) )
        factor = 1e15 ;
    elseif( str.equals( java.lang.String('T') ) )
        factor = 1e12 ;
    elseif( str.equals( java.lang.String('G') ) )
        factor = 1e9 ;
    elseif( str.equals( java.lang.String('M') ) )
        factor = 1e6 ;
    elseif( str.equals( java.lang.String('k') ) )
        factor = 1e3 ;
    elseif( str.equals( java.lang.String('h') ) )
        factor = 1e2 ;
    elseif( str.equals( java.lang.String('da') ) )
        factor = 1e1 ;
    elseif( str.equals( 'none' ) || str.equals( java.lang.String( '' ) ) )
        factor = 1e0 ;
    elseif( str.equals( java.lang.String('d') ) )
        factor = 1e-1 ;
    elseif( str.equals( java.lang.String('c') ) )
        factor = 1e-2 ;
    elseif( str.equals( java.lang.String('m') ) )
        factor = 1e-3 ;
    elseif( str.equals( java.lang.String( java.lang.Character.toChars( hex2dec( '03BC' ) ) ) ) )
        factor = 1e-6 ;
    elseif( str.equals( java.lang.String('n') ) )
        factor = 1e-9 ;
    elseif( str.equals( java.lang.String('p') ) )
        factor = 1e-12 ;
    elseif( str.equals( java.lang.String('f') ) )
        factor = 1e-15 ;
    elseif( str.equals( java.lang.String('a') ) )
        factor = 1e-18 ;
    elseif( str.equals( java.lang.String('z') ) )
        factor = 1e-21 ;
    elseif( str.equals( java.lang.String('y') ) )
        factor = 1e-24 ;
    else
        factor = 1 ;
        warning( 'UNIT-MODIFIER didn''t match ENG notation!' ) ; %#ok<WNTAG>
    end ;
end























