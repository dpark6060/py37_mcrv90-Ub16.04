function plotOverlay(R,varargin)

% if numel(size(R)) < 6, error('Input matrix must be 6-D'); end

%% parameters
par.bgthres = [0 nan];      % nan means auto-determine
par.selZ = [];              % empty means all positions, integer means numbers, float is stage positions
par.bgwl = nan;             % nan means last wavelength
par.xspac = 10;
par.yspac = 10;
par.xfac = 1;       % aspect ratio of axes
par.yfac = 1;
par.cbwidth = 20;
par.fgthres = [0 nan];
par.fgathres = uint16([0 100]); % either inabsolute (double) or percent (integer)
par.fginvert = 0;
par.selCmp = 1;
par.fgcmap = jet;
par.cbmode = 3;
par.bgcbar = 0;
par.autoclose = 0;
par.savename = [];          % empty means no saving
par.saveformat = 'png';
par.showtext = 1;

% Unmixed array
U = [];
if numel(varargin) >= 1
    U = varargin{1};
end

datainfo = [];
if numel(varargin) >= 3
    datainfo = varargin{3};
    if (isempty(par.savename))
        par.savename = datainfo.Name;
    end
end

% Copy parameters from input struct
if numel(varargin) >= 2
    cpar = varargin{2};
    % transfer all fields to parameter array
    fx = fieldnames(cpar);
    for j = 1:numel(fx)
        par = setfield(par,fx{j},getfield(cpar,fx{j}));
    end
    clear cpar j fx;
end


%% check selectors
% Z Selector must be in integer numbers, not stage position
if (isempty(par.selZ)), par.selZ = uint16(1:size(R,4)); end;
if (isfloat(par.selZ)) 
    oZ = par.selZ; par.selZ = [];
    for j = 1:numel(oZ)
        par.selZ = [par.selZ find(oZ(j) == datainfo.ZPositions)]; 
    end
    clear j oZ;
end

% bg selector
if (isnan(par.bgwl)), par.bgwl = numel(datainfo.Wavelengths); end;
if (max(par.bgwl) > 100) 
    par.bgwl = find(par.bgwl == datainfo.Wavelengths); 
end


%% filter inputs
R = R(:,:,:,par.selZ,:,par.bgwl,:);
if (~isempty(U)) U = U(:,:,:,par.selZ,:,par.selCmp); end
par.zpos = datainfo.ZPositions(par.selZ);
par.wls = datainfo.Wavelengths(par.bgwl);

if ( ~isempty(U) && ((numel(size(R)) ~= numel(size(U))) || any(size(U) ~= size(R))) )
    error('Size mismatch between foreground and background');
end

%% thresholds
if (isnan(par.bgthres(1))) par.bgthres(1) = min(R(:)); end
if (isnan(par.bgthres(2))) par.bgthres(2) = max(R(:)); end
if ~isempty(U) 
    if par.fginvert, U = -U; end;
    if (isnan(par.fgthres(1))) par.fgthres(1) = min(U(:)); end
    if (isnan(par.fgthres(2))) par.fgthres(2) = max(U(:)), end 
    if (~isfloat(par.fgthres(1))) 
        par.fgthres = double(par.fgthres);
        par.fgthres(1) = diff([0 max(U(:))])*par.fgthres(1)/100; 
        par.fgthres(2) = diff([0 max(U(:))])*par.fgthres(2)/100; 
    end
end

%% plot
for sl = 1:size(R,4)
    f=figure;
    whitebg(f);
    ps = get(f,'Position');
    ps(3) = par.xspac*(2.5+par.bgcbar*0.5)+(size(R,2)*par.xfac)+par.cbwidth*par.bgcbar+par.cbwidth;
    ps(4) = 2*par.yspac+round(size(R,1)*par.yfac);
    set(f,'Position',ps);
%     set(f,'Renderer','OpenGL');
%     set(f,'RendererMode','auto');
%     get(f);
    
    ax = axes('Units','pixel','Position',[par.xspac par.yspac round(size(R,2)*par.xfac) round(size(R,1)*par.yfac)]);
    
    % bg image
    img = R(:,:,1,sl,1);
    img = img - par.bgthres(1); img(img < 0) = 0;
    img = img ./ (par.bgthres(2) - par.bgthres(1)); img(img > 1) = 1;
    img_c = ind2rgb(round(img*64),gray);
    image(img_c);
    axis off;
    % bg colorbar
    if par.bgcbar
        par.cbarlabel{1} = num2str(par.bgthres(1),'%i');
        par.cbarlabel{2} = num2str(par.bgthres(2),'%i');
         plotColorbar(...
            [par.xspac*1.5+size(R,2),par.yspac,par.cbwidth,size(R,1)],...
            gray,par);
    end
    
    
    % fg image
    if ~isempty(U)
        fgimg = U(:,:,1,sl,1);
        fgimg = fgimg - par.fgthres(1); fgimg(fgimg < 0) = 0;
        fgimg = fgimg ./ (par.fgthres(2)- par.fgthres(1)); fgimg(fgimg > 1) = 1;
        fgimg_c = ind2rgb(round(fgimg*64),par.fgcmap);
        axes(ax);
        hold on;
        fgh = image(fgimg_c);
        axis off;
        
        % alpha
        if ~isfloat(par.fgathres)
            par.fgathres = double(par.fgathres);
            par.fgathres(1) = par.fgthres(1)+(double(par.fgathres(1))/100*diff(par.fgthres));
            par.fgathres(2) = par.fgthres(1)+(double(par.fgathres(2))/100*diff(par.fgthres));
        end
        aimg = U(:,:,1,sl,1);
        aimg = aimg - par.fgathres(1); aimg(aimg < 0) = 0;
        aimg = aimg ./ (par.fgathres(2)- par.fgathres(1)); aimg(aimg > 1) = 1;
        set(fgh,'AlphaData',aimg);
        
        astart = round((par.fgathres(1)-par.fgthres(1))/diff(par.fgthres)*64);
        if astart < 1, astart = 1; end;
        aend = round((par.fgathres(2)-par.fgthres(1)/diff(par.fgthres)*64));
        if(aend > 64), aend = 64; end;
        par.fgcmap(:,4) = ones(size(par.fgcmap,1),1);
        par.fgcmap(aend:astart,4) = linspace(1,0,astart-aend+1);
        par.fgcmap(astart:end,4) = 0;
        par.cbarlabel{1} = num2str(round(par.fgthres(1)),'%i');
        par.cbarlabel{2} = num2str(round(par.fgthres(2)),'%i');
        plotColorbar(...
            [par.xspac*(1.5+par.bgcbar*0.5)+round(size(R,2)*par.xfac)+par.cbwidth*par.bgcbar,par.yspac,par.cbwidth,round(size(R,1)*par.yfac)],...
            par.fgcmap,par);
    end
    
    axes(ax);
    if par.showtext
        text(4,size(R,2)-10,sprintf(' z: %.1fmm',par.zpos(sl)),'Backgroundcolor',[0 0 0],'Color',[1 1 1],'Margin',4);
    end
    
    if (isstruct(datainfo) && par.showtext)
        text(4,4,[' ' datainfo.Name],'Backgroundcolor',[0 0 0],'Color',[1 1 1],'Margin',4,'FontSize',10,'Interpreter','none','VerticalAlign','top');
    end
    
    % save
    if (~isempty(par.savename))
        pause(1);
        fr = getframe(f);
        imwrite(fr.cdata,[par.savename '_cmp' num2str(par.selCmp) '_z' strrep(num2str(par.zpos(sl),'%05.1f'),'.','-') 'mm.' par.saveformat]);
        
%         figure;image(fr.cdata);
    end
    
    % close if requested
    if par.autoclose,        close(f);    end
    
end