function [hdl,cfg] = nplot_error(varargin)
%
% NPLOT_ERROR plots data with error bars in various formats and offers the
% choice betweeen within-subject and between-subject errorbars
%
% USAGE: [hdl,cfg] = nplot_error('option','value',...)
%
% nplot_error take option value pairs as inputs and returns the axis handle
% for further manual plot commands and the config struct (cfg), which can
% be also edited manually and used as input for the next nplot_error call.
%
% INPUTS
% OPTION           (POSSIBLE) VALUES
% 'config'         struct with config options (from previous nplot_error)
%
% 'type'           'bar','dot','line','fact'
%
% 'data'           nSubjects*nConditions data matrix
% 'mean'           vector of (manually) computed means
% 'error'          vector of (manually) computed errors
%
% 'errorclass'     'bw' (between-subject) or 'ws' (within-subject)
% 'errortype'      'sem','sd','ci' (90%), 'none'
% 'color'          cell array of color specs (1*nCondition) [bar,dot,line] OR
%                  vector of dimension 1*nCondition with values between
%                  1:8 referring to the spcialized color palette defined in
%                  cfg.colorder (1=blue, 2=orange, 3=red, 4=green, 5=violet,
%                  6=yellow, 7=magenta, 8=brown)
%                  Colors can be specified using standard MATLAB characters
%                  (rgbcmykw) or as RGB vectors with values between 0 and 1
% 'errorcolor'     string or numerical color spec for errorbars
%                  special string 's' for "same as bar/dot"
%
% 'indivdata'      plot individual data points as little dots ('on','off')
% 'subjectnum'     plot subject numbers next to dots for individual data
%                  ('on','off')
%
% 'group'          cell array of vectors of subject numbers (rows in data
%                  matrix), specify xpos, xtick, and xlabel for one group
%                  only, it will be automatically replicated for the remaining
%                  groups
% 'gnames'         cell array of names for groups
%
% 'factor'         vector of number of level per factor, left-most factor
%                  rotates slowest [type=fact,prod(factor)==size(data,2)]
% 'factornames'    cell array of factor names [type=fact]
% 'levelnames'     cell array of cell array of levelnames (for each factor)
%                  [type=fact]
%
% 'xpos'           positions on x-axis for plotting the data, use this for
%                  visual grouping
% 'xtick'          ticks on x-axis for adding labels for data
% 'xlabel'         cell array of labels for xtick [type=bar/dot/line]
%                  (must correspond to length(xtick)
% 'xaxislabel'     label for x-axis (type=bar/dot/line]
% 'yaxislabel'     label for y-axis
% 'title'          plot title
%
% 'legend'         flag to turn on legend ('on','off')
% 'legnames'       cell array for different colors [type=bar/dot/line]
% 'leglocation'    location for legend (e.g. 'nw','ne','se','sw')
%
% 'figcmd'         command for figure/axes/subplots
% _________________________________________________________________________
% (c) Jan Glaescher (glaescher@uke.de)

% Hard-wired defaults in "cfg.defs" config struct (can be changed in the
% returned "cfg" struct and then re-submitted using
% nplot_error('config',cfg)
%
% cfg.defs.marker - marker for dot/line plot ('o')
% cfg.defs.markersize - marker size for dot/line plot default (10)
% cfg.defs.markeredgecolor - border color of marker ('k')
% cfg.defs.markerfacecolor - face color of marker ('w')
% cfg.defs.linewidth - line width for dot/line plot (1)
% cfg.defs.linestyle - linestyle for dot/line plot ('-')
% cfg.defs.fontsize - fontsize for ticklabels (14)
% cfg.defs.axisfontsize - fontsize for axis labels (16)
% cfg.defs.titlefontsize  fontsize for title (18)
% cfg.defs.colorder - order of colors for auto-assignment in type 'fact'
% cfg.defs.symorder - order of marker symbols for type 'fact'
% cfg.defs.lineorder - order of linestyles for type 'fact'
% cfg.defs.ebarhtick - horiz. tick marks at the end of the errorbar (0/1)
% _________________________________________________________________________
% nplot_error offers a new plot type: 'fact'orial, which parcellates a data
% matrix according to a factor level specification ('factor'). The
% left-most factor rotates slowest. The right-most factor rotates fastest and
% is plotted on the x-axis. Other factors are shown as different color or
% symbols and linestyles, which are automatically assigned. This factor
% specification must be valid for the data matrix specified in 'data'.
% Factorial plots automatically use a within-subject error. Ticklabels,
% axis labels, and legend labels are automatically extracted from the
% specifications in 'factornames' and 'levelnames'.
%
% If a within-subject error is also needed for one of the other plot types,
% then it needs to be explicitly selected ('errorclass','ws'). In addition,
% the 'factor' option has to be specified by usually supplying a single
% number ['factor',size(data,2)].
%
% 'xpos' can be used for visual grouping along the x-axis, e.g.
% 'xpos',[1:3 5:7] leaves a small gap between the first 3 and the second 3
% conditions. If the option 'group' is specified, then entire plot is
% "replicated" for each group along the x-axis. The chosen errorclass
% (between- or within-subjects) is used for each of the groups.
% _________________________________________________________________________
% Within-subject errors are computed according to Morey (2008) [3]. Data
% are first normalized w.r.t. to the subject's mean and then standard SEMs
% are computed on the normalized data. These SEM/SD/CI are corrected for
% the number of factor levels according to Morey (2008). While these
% errorbars are not identical to Loftus & Masson (2004) [1], as detailed in
% Franz & Loftus (2012) [3], they do provide a good approximation, which allows
% for different errorbars for each condition (the Loftus & Masson errorbars
% use the pooled error term and hence errorbar are identical for all
% levels). The method of paired differences (Franz & Loftus, 2012) [4]
% leaves some ambiguities, which errorbar should be used for a condition as
% all conditions occur in several paired differences. For a full discussion
% see
% [1] Loftus & Masson (1994). Using confidence intervals in within-subject
%     designs, Psychon Bull Rev, 1(4), 476-490.
% [2] Cousineau (2005), Confidence intervals in within-subject designs: a
%     simpler solution to Loftus & Masson's method, Tutorials in Quantitative
%     Psychology, 1(10, 42-45.
% [3] Morey (2008). Intervals from Normalized Data: A correction to
%     Cousineau (2005), Tutorials in Quantitative Psychology, 4(2), 61-64
% [4] Franz & Loftus (2012). Standard errors and confidence intervals in
%     within-subject designs: Generalizing Loftus & Masson (1994) and avoiding
%     the biases of alternative accounts, Psychon Bull Rev, 19, 395-404

%% preliminaries
% check for inverse cumulative T-distribution (for computing CIs)
has_t(1) = ~isempty(which('spm_invTcdf')); % from SPM distribution
has_t(2) = ~isempty(which('tinv'));

if mod(nargin,2) ~= 0
  error('Inconsistent number of option-value pairs')
end

%% parse input
inargs = varargin;
cfg = parse_input(inargs{:});

%% convert data matrix in means and errors (depending on type)
if ~isempty(cfg.data)
  if isfield(cfg,'group')
    for g = 1:length(cfg.group)
      cfg.mean(g,:) = nanmean(cfg.data(cfg.group{g},:));
    end
  else
    cfg.mean = nanmean(cfg.data);
  end
  if strcmp(cfg.errorclass,'bw')
    if strcmp(cfg.errortype,'sem')
      if isfield(cfg,'group')
        for g = 1:length(cfg.group)
          cfg.error(g,:) = nansem(cfg.data(cfg.group{g},:));
        end
      else
        cfg.error = nansem(cfg.data);
      end
    elseif strcmp(cfg.errortype,'sd')
      if isfield(cfg,'group')
        for g = 1:length(cfg.group)
          cfg.error(g,:) = nanstd(cfg.data(cfg.group{g},:));
        end
      else
        cfg.error = nanstd(cfg.data);
      end
    elseif strcmp(cfg.errortype,'ci')
      j = find(has_t,1);
      if ~isempty(j)
        if j == 1
          if isfield(cfg,'group')
            for g = 1:length(cfg.group)
              m(g) = spm_invTcdf(0.95,size(cfg.data(cfg.group{1},:),1)-1);
            end
          else
            m = spm_invTcdf(0.95,size(cfg.data,1)-1);
          end
        elseif j == 2
          if isfield(cfg,'group')
            for g = 1:length(cfg.group)
              m(g) = tinv(0.95,size(cfg.data(cfg.group{1},:),1)-1);
            end
          else
            m = tinv(0.95,size(cfg.data,1)-1);
          end
        end
        if isfield(cfg,'group')
          for g = 1:length(cfg.group)
            cfg.error(g,:) = nanstd(cfg.data(cfg.group{g},:)) * m(g);
          end
        else
          cfg.error = nanstd(cfg.data) * m;
        end
      else
        fprintf(2,'Warning: CIs cannot be computed (missing spm_invTcdf ot tinv). Using SEM instead.');
        cfg.error = nansem(cfg.data);
      end
    end
  elseif strcmp(cfg.errorclass,'ws')
    % normalize data
    if isfield(cfg,'factor')
      nLevels = prod(cfg.factor);
    else
      nLevels = size(cfg.data,2);
    end
    
    if isfield(cfg,'group')
      for g = 1:length(cfg.group)
        ndata = cfg.data(cfg.group{g},:) - ( repmat(nanmean(cfg.data(cfg.group{g},:),2),1,...
          size(cfg.data(cfg.group{g},:),2)) - repmat(nanmean(cfg.data(cfg.group{g},:)),...
          size(cfg.data(cfg.group{g},:),1),1) );
        if strcmp(cfg.errortype,'sem')
          cfg.error(g,:) = nansem(ndata) * sqrt( nLevels / (nLevels-1) );
        elseif strcmp(cfg.errortype,'sd')
          cfg.error(g,:) = nanstd(ndata) * sqrt( nLevels / (nLevels-1) );
        elseif strcmp(cfg.errortype,'ci')
          j = find(has_t,1);
          if ~isempty(j)
            if j == 1
              m(g) = spm_invTcdf(0.95,size(cfg.data(cfg.group{g},:),1)-1);
            elseif j == 2
              m(g) = tinv(0.95,size(cfg.data(cfg.group{g},:),1)-1)
            end
            cfg.error(g,:) = nanstd(ndata) * m(g) * sqrt( nLevels / (nLevels-1) );
          else
            fprintf(2,'Warning: CIs cannot be computed (missing spm_invTcdf ot tinv). Using SEM instead.');
            cfg.error(g,:) = nansem(ndata) * sqrt( nLevels / (nLevels-1) );
          end
        end
        
      end
    else
      ndata = cfg.data - ( repmat(nanmean(cfg.data,2),1,size(cfg.data,2)) - ...
        repmat(nanmean(cfg.data(:)),size(cfg.data)) );
      if strcmp(cfg.errortype,'sem')
        cfg.error = nansem(ndata) * sqrt( nLevels / (nLevels-1) );
      elseif strcmp(cfg.errortype,'sd')
        cfg.error = nanstd(ndata) * sqrt( nLevels / (nLevels-1) );
      elseif strcmp(cfg.errortype,'ci')
        j = find(has_t,1);
        if ~isempty(j)
          if     j == 1; m = spm_invTcdf(0.95,size(cfg.data,1)-1);
          elseif j == 2; m = tinv(0.95,size(cfg.data,1)-1); end
          cfg.error = nanstd(ndata) * m * sqrt( nLevels / (nLevels-1) );
        else
          fprintf(2,'Warning: CIs cannot be computed (missing spm_invTcdf ot tinv). Using SEM instead.');
          cfg.error = nansem(cfg.data);
        end
      end
    end
  elseif strcmp(cfg.errorclass,'none')
    cfg.error = zeros(size(cfg.mean));
  end
end

%% === Start plotting
% --- determine figure window
if isfield(cfg,'figcmd')
  fig = eval(cfg.figcmd);
  hdl = gca;
else
  fig = figure;
  hdl = axes;
  set(hdl,'fontsize',cfg.defs.fontsize);
end
hold on;

if isfield(cfg,'group')
  nGroup = length(cfg.group);
  for g = 2:nGroup
    cfg.xpos(g,:) = cfg.xpos(g-1,:) + cfg.xpos(g-1,end) + 3;
    cfg.xtick(g,:) = cfg.xtick(g-1,:) + cfg.xpos(g-1,end) + 3;
    cfg.xlabel(g,:) = cfg.xlabel(g-1,:);
  end
else
  nGroup = 1;
end

switch cfg.type
  
  % --- bar plot: (1) bars and (2) errorbars
  case 'bar'
    for x = 1:nGroup
      for j = 1:length(cfg.mean(x,:))
        b(x,j) = bar(cfg.xpos(x,j),cfg.mean(x,j));
        b(x,j).FaceColor = cfg.color{j};
        %set(b(x,j),'facecolor',cfg.color{j});
        e(x,j) = line([cfg.xpos(x,j) cfg.xpos(x,j)],...
          [cfg.mean(x,j)-cfg.error(x,j) cfg.mean(x,j)+cfg.error(x,j)]);
        e(x,j).Color     = cfg.errorcolor{j};
        e(x,j).LineWidth = cfg.defs.linewidth;
        %set(e(x,j),'color',cfg.errorcolor{j},'linewidth',cfg.defs.linewidth);
        if cfg.defs.ebarhtick
          len = max(length(cfg.xpos(x,:))/100,0.1);
          up(x,j) = line([cfg.xpos(x,j)-len cfg.xpos(x,j)+len],...
            [cfg.mean(x,j)+cfg.error(x,j) cfg.mean(x,j)+cfg.error(x,j)]);
          lo(x,j) = line([cfg.xpos(x,j)-len cfg.xpos(x,j)+len],...
            [cfg.mean(x,j)-cfg.error(x,j) cfg.mean(x,j)-cfg.error(x,j)]);
          up(x,j).Color     = cfg.errorcolor{j};
          up(x,j).LineWidth = cfg.defs.linewidth;
          lo(x,j).Color     = cfg.errorcolor{j};
          lo(x,j).LineWidth = cfg.defs.linewidth;
          %set(up(x,j),'color',cfg.errorcolor{j},'linewidth',cfg.defs.linewidth);
          %set(lo(x,j),'color',cfg.errorcolor{j},'linewidth',cfg.defs.linewidth);
        end
      end
    end
    b(1).BaseLine.Visible = 'off';
    yl = hdl.YLim;
    yl(1) = yl(1) - range(yl)/10;
    hdl.YLim = yl;
    
    
    % --- dot/line plot: (1) errorbars and (2) dots/lines
  case { 'dot', 'line' }
    for g = 1:nGroup
      if strcmp(cfg.type,'line')
        if ~isequalcell(cfg.color)
          fprintf(2,'Warning: Cannot plot line plot with different color specs. Using dot plot instead.\n');
        else
          li = plot(cfg.xpos(g,:),cfg.mean(g,:));
          li.Marker    = 'none';
          li.Color     = cfg.color{1};
          li.LineStyle = '-';
          li.LineWidth = cfg.defs.linewidth;
        end
      end
      for j = 1:length(cfg.mean(g,:))
        e(g,j) = line([cfg.xpos(g,j) cfg.xpos(g,j)],...
          [cfg.mean(g,j)-cfg.error(g,j) cfg.mean(g,j)+cfg.error(g,j)]);
        e(g,j).Color = cfg.errorcolor{j};
        e(g,j).LineWidth = cfg.defs.linewidth;
        %set(e(g,j),'color',cfg.errorcolor{j},'linewidth',cfg.defs.linewidth);
        if cfg.defs.ebarhtick
          len = max(length(cfg.xpos(g,:))/100,0.1);
          up(g,j) = line([cfg.xpos(g,j)-len cfg.xpos(g,j)+len],...
            [cfg.mean(g,j)+cfg.error(g,j) cfg.mean(g,j)+cfg.error(g,j)]);
          lo(g,j) = line([cfg.xpos(g,j)-len cfg.xpos(g,j)+len],...
            [cfg.mean(g,j)-cfg.error(g,j) cfg.mean(g,j)-cfg.error(g,j)]);
          up(g,j).Color     = cfg.errorcolor{j};
          up(g,j).LineWidth = cfg.defs.linewidth;
          lo(g,j).Color     = cfg.errorcolor{j};
          lo(g,j).LineWidth = cfg.defs.linewidth;
          %set(up(g,j),'color',cfg.errorcolor{j},'linewidth',cfg.defs.linewidth);
          %set(lo(g,j),'color',cfg.errorcolor{j},'linewidth',cfg.defs.linewidth);
        end
        if strcmp(cfg.defs.markerfacecolor,'s')
          facecol = cfg.color{j};
        else
          facecol = 'w';
        end
        b(j) = plot(cfg.xpos(g,j),cfg.mean(g,j));
        b(j).Marker = 'o';
        b(j).LineStyle = 'none';
        b(j).MarkerSize = cfg.defs.markersize;
        b(j).MarkerEdgeColor = cfg.color{j};
        b(j).MarkerFaceColor = facecol;
      end
    end
    
    % --- factorial plot with ws-errors (always as lin plot)
  case 'fact'
    for g=1:nGroup
      % --- sort data into factor levels
      % number of factor levels of the fastest rotating factor
      if mod(length(cfg.mean(g,:)) ./ length(cfg.xpos(g,:)),1)
        error('Number of effects is not a multiple of number of x-positions.')
      end
      nx = length(cfg.xpos(g,:));
      nl = length(cfg.mean(g,:)) ./ length(cfg.xpos(g,:));
      M = zeros(nl,nx,nGroup); E = zeros(nl,nx,nGroup);
      for f=1:nl
        M(f,:,g) = cfg.mean(g,((f-1)*nx)+(1:nx));
        E(f,:,g) = cfg.error(g,((f-1)*nx)+(1:nx));
      end
      % create lookup table for colors, linestyles, symbols
      lut = create_lut(cfg.factor);
      cx = strcmp('color',cfg.attrib);
      sx = strcmp('syms',cfg.attrib);
      if mod(size(M,1),2) % uneven number of xticks
        h = (size(M,1)-1) / 2;
        offx = (-h:h) * 0.05;
      else % even number of xticks
        h = size(M,1) / 2;
        offx = [ -h:-1 1:h ] * 0.05;
      end
      for j=1:size(M,1)
        lw = cfg.defs.linewidth;
        if any(sx)
          lord = cfg.defs.lineorder{lut(j,sx)};
        else
          lord = cfg.defs.lineorder{1};
        end
        li(g,j) = plot(cfg.xpos(g,:)+offx(j),M(j,:,g));
        li(g,j).Color = cfg.color{j};
        li(g,j).LineStyle = lord;
        li(g,j).LineWidth = cfg.defs.linewidth;
        li(g,j).Marker = 'none';
        for x = 1:size(M,2)
          e(g,j,x) = line([cfg.xpos(g,x)+offx(j) cfg.xpos(g,x)+offx(j)],...
            [M(j,x,g)-E(j,x,g) M(j,x,g)+E(j,x,g)]);
          e(g,j,x).Color = cfg.errorcolor{j};
          e(g,j,x).LineWidth = cfg.defs.linewidth;
          %set(e(g,j,x),'color',cfg.errorcolor{j},'linewidth',cfg.defs.linewidth);
          if cfg.defs.ebarhtick
            len = max(length(cfg.xpos(g,:))/100,0.05);
            up(g,j,x) = line([cfg.xpos(g,x)+offx(j)-len cfg.xpos(g,x)+offx(j)+len],...
              [M(j,x,g)+E(j,x,g) M(j,x,g)+E(j,x,g)]);
            lo(g,j,x) = line([cfg.xpos(g,x)+offx(j)-len cfg.xpos(g,x)+offx(j)+len],...
              [M(j,x,g)-E(j,x,g) M(j,x,g)-E(j,x,g)]);
            up(g,j,x).Color     = cfg.errorcolor{j};
            up(g,j,x).LineWidth = cfg.defs.linewidth;
            lo(g,j,x).Color     = cfg.errorcolor{j};
            lo(g,j,x).LineWidth = cfg.defs.linewidth;
            %set(up(g,j,x),'color',cfg.errorcolor{j},'linewidth',cfg.defs.linewidth);
            %set(lo(g,j,x),'color',cfg.errorcolor{j},'linewidth',cfg.defs.linewidth);
          end
          p(g,j,x) = plot(cfg.xpos(g,:)+offx(j),M(j,:,g));
          if any(sx); so = lut(j,sx); else so = 1; end
          if ~isempty(lut); co = lut(j,cx); else co = 1; end
          p(g,j,x).LineStyle = 'none';
          p(g,j,x).Marker = cfg.defs.symorder{so};
          p(g,j,x).MarkerSize = cfg.defs.markersize;
          p(g,j,x).MarkerEdgeColor = cfg.color{j};
          p(g,j,x).MarkerFaceColor = cfg.defs.markerfacecolor;
          %set(p(g,j,x),...
          %  'linestyle','none',...
          %  'marker',cfg.defs.symorder{so},...
          %  'markersize',cfg.defs.markersize,...
          %  'markeredgecolor', cfg.color{j},...
          %  'markerfacecolor', cfg.defs.markerfacecolor);
          %'markeredgecolor', cfg.defs.colorder{co},...
        end
      end
    end
end

%% --- individual data points and labels
if strcmp(cfg.indivdata,'on')
  if isfield(cfg,'data')
    for g = 1:nGroup
      for j = 1:size(cfg.data(g,:),2)
        for s = 1:size(cfg.data(cfg.group{g},:),1)
          if strcmp(cfg.subjectnum,'on')
            if mod(s,2)
              plot(cfg.xpos(g,j)-0.07,cfg.data(cfg.group{g}(s),j),'k.',...
                'markersize',cfg.defs.indivdatamarker)
              text(cfg.xpos(g,j)-0.1,cfg.data(cfg.group{g}(s),j),num2str(cfg.group{g}(s)),...
                'fontsize',8,'horizontalalignment','right')
            else
              plot(cfg.xpos(g,j)+0.07,cfg.data(cfg.group{g}(s),j),'k.',...
                'markersize',cfg.defs.indivdatamarker)
              text(cfg.xpos(g,j)+0.1,cfg.data(cfg.group{g}(s),j),num2str(cfg.group{g}(s)),...
                'fontsize',8,'horizontalalignment','left')
            end
          else
            plot(cfg.xpos(g,j)+0.07,cfg.data(cfg.group{g}(s),j),'k.',...
              'markersize',cfg.defs.indivdatamarker)
          end
        end
      end
    end
  else
    fprintf(2,'Cannot print individual data points. cfg.data is missing.\n');
  end
end

%% --- general axis settings

if nGroup > 1
  % compute xposition for group labels (used later)
  for g = 1:nGroup
    gx(g) = cfg.xpos(g,1) + range(cfg.xpos(g,:))/2;
  end
  xtick = cfg.xtick';
  xtick = xtick(:)';
  xpos = cfg.xpos';
  xpos = xpos(:)';
  xlabels = cfg.xlabel';
  xlabels = xlabels(:)';
  xlim = [(min(xpos)-1) (max(xpos)+1)];
else
  if isfield(cfg,'xtick')
    xtick = cfg.xtick;
  else
    xtick = cfg.xpos;
  end
  xpos = cfg.xpos;
  xlabels = cfg.xlabel;
  xlim = [min(xpos)-1 max(xpos)+1];
end

set(hdl,'xtick',xtick,'xticklabel',xlabels,...
  'xlim',xlim,'fontsize',cfg.defs.fontsize);
if isfield(cfg,'xaxislabel')
  xlabel(cfg.xaxislabel,'fontsize',cfg.defs.axisfontsize);
elseif strcmp(cfg.type,'fact')
  xlabel(cfg.factornames{end},'fontsize',cfg.defs.axisfontsize);
end
if isfield(cfg,'yaxislabel')
  ylabel(cfg.yaxislabel,'fontsize',cfg.defs.axisfontsize);
end
if isfield(cfg,'title')
  title(cfg.title,'fontsize',cfg.defs.titlefontsize);
end

%% --- legends
if strcmp(cfg.legend,'on')
  if strcmp(cfg.type,'fact')
    idx = strcmp('color',cfg.attrib);
    legend(li,cfg.levelnames{idx},'location',cfg.leglocation)
    legend boxoff
  elseif any(strcmp(cfg.type,{'bar','dot','line'}))
    if iscellstr(cfg.color)
      u = unique(cfg.color,'stable');
      for f=1:length(u)
        idx(f) = find(strcmp(u{f},cfg.color),1,'first');
      end
    else
      u = unique(cat(1,cfg.color{:}),'rows','stable');
      for f=1:size(u,1)
        for g = 1:length(cfg.color)
          if u(f,:) == cfg.color{g}
            idx(f) = g;
            break;
          end
        end
      end
    end
    if length(idx) == length(cfg.legnames)
      legend(b(1,idx),cfg.legnames,'location',cfg.leglocation);
      legend boxoff
    else
      error('Mismatch between unique colors and legend names')
    end
  end
end

%% plot group labels if necessary
if nGroup > 1
  yl = get(gca,'ylim');
  ypos = yl(2) - (yl(2)-yl(1))/50;
  for g=1:nGroup
    if isfield(cfg,'gnames')
      tx(g) = text(gx(g),ypos,cfg.gnames{g});
      tx(g).FontSize = cfg.defs.fontsize;
      tx(g).FontWeight = 'bold';
      tx(g).HorizontalAlignment = 'c';
      tx(g).VerticalAlignment = 't';
      %  'fontweight','bold','horizontalalignment','c','verticalalignment','t');
    else
      tx(g) = text(gx(g),ypos,sprintf('Group %d',g));
      tx(g).Fontsize = cfg.defs.fontsize;
      tx(g).FontWeight = 'bold';
      tx(g).HorizontalAlignment = 'c';
      tx(g).VerticalAlignment = 't';
      %  'fontweight','bold','horizontalalignment','c','verticalalignment','t');
    end
  end
end
return

%% nested sub-functions
  function cfg = parse_input(varargin)
    
    if strcmp(varargin{1},'config')
      cfg = varargin{2};
      return;
    end
    
    cfg.defs.marker = 'o';
    cfg.defs.markersize = 8;
    cfg.defs.markeredgecolor = 'k';
    cfg.defs.markerfacecolor = 'w';
    cfg.defs.linewidth = 2;
    cfg.defs.linestyle = '-';
    cfg.defs.fontsize = 14;
    cfg.defs.axisfontsize = 14;
    cfg.defs.titlefontsize = 16;
    cfg.defs.colorder = {[0.1211    0.4648    0.7031], ...
      [0.9961    0.4961    0.0586], ...
      [0.8359    0.1523    0.1562], ...
      [0.1719    0.6250    0.1719], ...
      [0.5781    0.4023    0.7383], ...
      [1.0000    0.7500         0], ...
      [0.8945    0.2578    0.4570], ...
      [0.4141    0.2031    0.1680]};
    
    cfg.defs.symorder = {'o','^','d','x'};
    cfg.defs.lineorder = {'-','--','-.',':'};
    cfg.defs.ebarhtick = 1;
    
    cfg.defs.indivdatamarker = 8;
    
    i = 1;
    
    while i <= nargin-1
      
      arg = varargin{i};
      val = varargin{i+1};
      
      switch lower(arg)
        
        case 'config'
          
          cfg = val; % copy input config to output
          
        case 'type'
          
          if ~ischar(val)
            error('Value for option "%s" must be a character string',arg)
          end
          
          opt = {'bar','dot','line','fact'};
          o = strncmpi(val,opt,3);
          
          if ~any(o)
            error('Unknown value to option "%s".',arg)
          end
          
        case 'errorclass'
          
          if ~ischar(val)
            error('Value for option "%s" must be a character string',arg)
          end
          
          opt = {'bw','ws','none'};
          o = strncmpi(val,opt,2);
          
          if ~any(o)
            error('Unknown value to option "%s".',arg)
          end
          
        case 'errortype'
          
          if ~ischar(val)
            error('Value for option "%s" must be a character string',arg)
          end
          
          opt = {'sd','sem','ci'};
          o = strncmpi(val,opt,2);
          
          if ~any(o)
            error('Unknown value to option "%s".',arg)
          end
          
        case 'data'
          
          if ~isnumeric(val) || any(size(val) == 1)
            error('Value for option "%s" must be a matrix.',arg)
          end
          
        case { 'mean', 'error', 'factor', 'xtick', 'xpos' }
          
          if ~isnumeric(val) || ~any(size(val)) == 1
            error('Value for option "%s" must be a vector.',arg)
          end
          
        case { 'factornames', 'levelnames', 'group', 'xlabel', 'gnames' }
          
          if ~iscell(val)
            error('Value for option "%s" must be a cell array.',arg)
          end
          
        case { 'errorcolor', 'legnames', 'xaxislabel', 'yaxislabel', 'title', 'leglocation' }
          
          if isnumeric(val)
            error('Option "%s" requires either a character string or a cell array.',arg)
          end
          
        case 'color'
          
          if ~isvector(val) && ~iscell(val)
            error('Option "%s" requires either a vector or a cell array.',arg);
          end
          
        case { 'indivdata', 'subjectnum', 'legend' }
          
          if ~ ( strcmpi(val,'on') || strcmpi(val,'off') )
            error('Option "%s" requires ''on'' or ''off''.',arg)
          end
          
        case 'figcmd'
          
          if ~ischar(val)
            error('Option "%s" requires an character string (e.g. figure(2) or subplot(2,2,1)).',arg)
          end
          
        otherwise
          error('%s: unknown option string.',arg)
      end
      
      cfg.(lower(arg)) = val;
      i = i + 2;
      
    end
    
    % === some more validity checks ===
    % data/mean/error
    if isfield(cfg,'data')
      cfg.mean = [];
      cfg.error = [];
      nData = size(cfg.data,2);
    else
      cfg.data = [];
      if isfield(cfg,'mean')
        nData = length(cfg.mean);
        if isfield(cfg,'error')
          if length(cfg.mean) ~= length(cfg.error)
            error('Uequal lengths in "mean" and "error" vectors.')
          end
        else
          error('"Error" vector is missing.')
        end
      else
        error('"Mean" vector is missing. No data.')
      end
    end
    
    
    % setting for different plot types
    switch cfg.type
      
      case { 'bar', 'dot', 'line' }
        
        if ~isfield(cfg,'legend')
          cfg.legend = 0;
        end
        if ~isfield(cfg,'leglocation')
          cfg.leglocation = 'ne';
        end
        if ~isfield(cfg,'errorclass')
          cfg.errorclass = 'bw';
        end
        if ~isfield(cfg,'errortype')
          cfg.errortype = 'sem';
        end
        if isfield(cfg,'xpos')
          if length(cfg.xpos) ~= nData
            error('Mismatching number of x positions.')
          end
        else
          cfg.xpos = 1:nData;
        end
        if isfield(cfg,'color')
          if iscell(cfg.color)
            if length(cfg.color) > 1
              if length(cfg.color) ~= nData
                error('Mismatching number of color spects.')
              end
            else
              cfg.color = repmat(cfg.color,1,nData);
            end
          else % it's a vector
            if length(cfg.color) == 3 && any(mod(cfg.color,1))
              % it's a 1x3 color spec
              cfg.color = repmat({cfg.color},1,nData);
            else
              if length(cfg.color) ~= nData
                error('Mismatching number of color spects.')
              end
              if all(mod(cfg.color,1)==0) && all(cfg.color <= 8)
                % it's a vector of integers -> use cfg.colorder
                tmp = cfg.color;
                cfg.color = cell(1,length(tmp));
                for f = 1:length(tmp)
                  cfg.color{f} = cfg.defs.colorder{tmp(f)};
                end
              else
                error('Incompatible color spec format.')
              end
            end
          end
        else
          if strcmp(cfg.type,'bar')
            cfg.color = repmat({[.9 .9 .9]},1,nData);
          else
            cfg.color = repmat(cfg.defs.colorder(1),1,nData);
          end
        end
        if isfield(cfg,'errorcolor')
          if length(cfg.errorcolor) == 1
            if ~iscell(cfg.errorcolor)
              if strcmp(cfg.errorcolor,'s')
                cfg.errorcolor = cfg.color;
              else
                cfg.errorcolor = repmat({cfg.errorcolor},1,nData);
              end
            else
              if strcmp(cfg.errorcolor{1},'s')
                cfg.errorcolor = cfg.color;
              else
                cfg.errorcolor = repmat(cfg.errorcolor,1,nData);
              end
            end
          else
            if length(cfg.errorcolor) ~= length(cfg.color)
              error('Mismatch between color specs for errorbars and means.');
            end
          end
        else
          if strcmp(cfg.type,'bar')
            cfg.errorcolor = repmat({'k'},1,nData);
          else
            cfg.errorcolor = cfg.color;
          end
        end
        if isfield(cfg,'xlabel')
          if isfield(cfg,'xtick')
            if length(cfg.xlabel) ~= length(cfg.xtick)
              error('Mismatch in number of xTicks and xTickLabels.')
            end
          else
            cfg.xtick = 1:length(cfg.xlabel);
          end
        else
          for l = 1:nData; cfg.xlabel{l} = num2str(l); end
        end
        if ~isfield(cfg,'indivdata')
          cfg.indivdata = 0;
        end
        if ~isfield(cfg,'subjectnum')
          cfg.subjectnum = 0;
        end
        
      case 'fact'
        
        if ~isfield(cfg,'legend')
          cfg.legend = 1;
        end
        if ~isfield(cfg,'leglocation')
          cfg.leglocation = 'ne';
        end
        if ~isfield(cfg,'errorclass')
          cfg.errorclass = 'ws';
        end
        if ~isfield(cfg,'errortype')
          cfg.errortype = 'sem';
        end
        if ~isfield(cfg,'factor')
          error('Missing required option "factor".')
        else
          if prod(cfg.factor) ~= nData
            error('Mismatch in factor specifications.')
          end
          if length(cfg.factor) == 1
            cfg.attrib = {'xpos'};
          elseif length(cfg.factor) == 2
            cfg.attrib = {'color','xpos'};
          elseif length(cfg.factor) == 3
            cfg.attrib = {'color','syms','xpos'};
          else
            error('Too many factors for visual features.')
          end
        end
        if isfield(cfg,'factornames')
          if length(cfg.factornames) ~= length(cfg.factor)
            error('Mismatch in factor names.')
          end
        else
          for t=1:length(cfg.factor)
            cfg.factornames{t} = char(64+t); % default names {'A','B','C',...}
          end
        end
        if isfield(cfg,'levelnames')
          for l=1:length(cfg.levelnames)
            if length(cfg.levelnames{l}) ~= cfg.factor(l)
              error('Mismatch in levelNames for factor %d',l)
            end
          end
        else
          for l = 1:length(cfg.levelnames)
            for g = 1:length(cfg.factor{l})
              cfg.levelnames{l}{g} = num2str(g);
            end
          end
        end
        if isfield(cfg,'color')
          idx = strcmp('color',cfg.attrib);
          if any(idx) && length(cfg.color) ~= cfg.factor(idx)
            error('Mismatch in color specs for factor.')
          end
        else
          idx = strcmp('color',cfg.attrib);
          if any(idx)
            for c=1:cfg.factor(idx)
              cfg.color{c} = cfg.defs.colorder{c};
            end
          else
            cfg.color = repmat(cfg.defs.colorder(1),1,cfg.factor(1));
          end
        end
        if isfield(cfg,'errorcolor')
          if length(cfg.errorcolor) == 1
            if ~iscell(cfg.errorcolor)
              if strcmp(cfg.errorcolor,'s')
                cfg.errorcolor = cfg.color;
              else
                cfg.errorcolor = repmat({cfg.errorcolor},1,length(cfg.color));
              end
            else
              if strcmp(cfg.errorcolor{1},'s')
                cfg.errorcolor = cfg.color;
              else
                cfg.errorcolor = repmat(cfg.errorcolor,1,length(cfg.color));
              end
            end
          else
            if length(cfg.errorcolor) ~= length(cfg.color)
              error('Mismatch between number of color specs for errorbars and means.')
            end
          end
        else
          cfg.errorcolor = cfg.color;
        end
        if isfield(cfg,'xpos')
          idx = strcmp('xtick',cfg.attrib);
          if length(cfg.xpos) ~= cfg.factor(idx)
            error('Mismatch in number of specified xpos for factor')
          end
        else
          idx = strcmp('xpos',cfg.attrib);
          cfg.xpos = 1:cfg.factor(idx);
        end
        if isfield(cfg,'xlabel')
          if length(cfg.xlabel) ~= length(cfg.xpos)
            error('Mismatch in specfied xTickLabels for factor.')
          end
        else
          idx = strcmp('xpos',cfg.attrib);
          cfg.xlabel = cfg.levelnames{idx};
        end
        if isfield(cfg,'gnames')
          if ~isfield(cfg,'group')
            error('Group names without group specification.')
          else
            if length(cfg.gnames) ~= length(cfg.group)
              error('Mismatch in number of group vectors and group names.')
            end
          end
          
        end
        if ~isfield(cfg,'indivdata')
          cfg.indivdata = 0;
        end
        if ~isfield(cfg,'subjectnum')
          cfg.subjectnum = 0;
        end
        
    end
    
  end

  function y = nanmean(x,dim)
    % FORMAT: Y = NANMEAN(X,DIM)
    %
    %    Average or mean value ignoring NaNs
    %
    %    This function enhances the functionality of NANMEAN as distributed in
    %    the MATLAB Statistics Toolbox and is meant as a replacement (hence the
    %    identical name).
    %
    %    NANMEAN(X,DIM) calculates the mean along any dimension of the N-D
    %    array X ignoring NaNs.  If DIM is omitted NANMEAN averages along the
    %    first non-singleton dimension of X.
    %
    %    Similar replacements exist for NANSTD, NANMEDIAN, NANMIN, NANMAX, and
    %    NANSUM which are all part of the NaN-suite.
    %
    %    See also MEAN
    
    % -------------------------------------------------------------------------
    %    author:      Jan Gl?scher
    %    affiliation: Neuroimage Nord, University of Hamburg, Germany
    %    email:       glaescher@uke.uni-hamburg.de
    %
    %    $Revision: 1.1 $ $Date: 2004/07/15 22:42:13 $
    
    if isempty(x)
      y = NaN;
      return
    end
    
    if nargin < 2
      dim = find(size(x)~=1, 1 );
      if isempty(dim)
        dim = 1;
      end
    end
    
    % Replace NaNs with zeros.
    nans = isnan(x);
    x(isnan(x)) = 0;
    
    % denominator
    count = size(x,dim) - sum(nans,dim);
    
    % Protect against a  all NaNs in one dimension
    i = find(count==0);
    count(i) = ones(size(i));
    
    y = sum(x,dim)./count;
    y(i) = i + NaN;
    
    
    
    % $Id: nanmean.m,v 1.1 2004/07/15 22:42:13 glaescher Exp glaescher $
    
    
  end

  function y = nanstd(x,dim,flag)
    % FORMAT: Y = NANSTD(X,DIM,FLAG)
    %
    %    Standard deviation ignoring NaNs
    %
    %    This function enhances the functionality of NANSTD as distributed in
    %    the MATLAB Statistics Toolbox and is meant as a replacement (hence the
    %    identical name).
    %
    %    NANSTD(X,DIM) calculates the standard deviation along any dimension of
    %    the N-D array X ignoring NaNs.
    %
    %    NANSTD(X,DIM,0) normalizes by (N-1) where N is SIZE(X,DIM).  This make
    %    NANSTD(X,DIM).^2 the best unbiased estimate of the variance if X is
    %    a sample of a normal distribution. If omitted FLAG is set to zero.
    %
    %    NANSTD(X,DIM,1) normalizes by N and produces the square root of the
    %    second moment of the sample about the mean.
    %
    %    If DIM is omitted NANSTD calculates the standard deviation along first
    %    non-singleton dimension of X.
    %
    %    Similar replacements exist for NANMEAN, NANMEDIAN, NANMIN, NANMAX, and
    %    NANSUM which are all part of the NaN-suite.
    %
    %    See also STD
    
    % -------------------------------------------------------------------------
    %    author:      Jan Gl?scher
    %    affiliation: Neuroimage Nord, University of Hamburg, Germany
    %    email:       glaescher@uke.uni-hamburg.de
    %
    %    $Revision: 1.1 $ $Date: 2004/07/15 22:42:15 $
    
    if isempty(x)
      y = NaN;
      return
    end
    
    if nargin < 3
      flag = 0;
    end
    
    if nargin < 2
      dim = find(size(x)~=1, 1 );
      if isempty(dim)
        dim = 1;
      end
    end
    
    
    % Find NaNs in x and nanmean(x)
    nans = isnan(x);
    avg = nanmean(x,dim);
    
    % create array indicating number of element
    % of x in dimension DIM (needed for subtraction of mean)
    tile = ones(1,max(ndims(x),dim));
    tile(dim) = size(x,dim);
    
    % remove mean
    x = x - repmat(avg,tile);
    
    count = size(x,dim) - sum(nans,dim);
    
    % Replace NaNs with zeros.
    x(isnan(x)) = 0;
    
    
    % Protect against a  all NaNs in one dimension
    i = find(count==0);
    
    if flag == 0
      y = sqrt(sum(x.*x,dim)./max(count-1,1));
    else
      y = sqrt(sum(x.*x,dim)./max(count,1));
    end
    y(i) = i + NaN;
    
    % $Id: nanstd.m,v 1.1 2004/07/15 22:42:15 glaescher Exp glaescher $
  end

  function y = nansem(x,dim)
    % FORMAT: Y = NANSEM(X,DIM)
    %
    %    Standard error of the mean ignoring NaNs
    %
    %    NANSTD(X,DIM) calculates the standard error of the mean along any
    %    dimension of the N-D array X ignoring NaNs.
    %
    %    If DIM is omitted NANSTD calculates the standard deviation along first
    %    non-singleton dimension of X.
    %
    %    Similar functions exist: NANMEAN, NANSTD, NANMEDIAN, NANMIN, NANMAX, and
    %    NANSUM which are all part of the NaN-suite.
    
    % -------------------------------------------------------------------------
    %    author:      Jan Gl?scher
    %    affiliation: Neuroimage Nord, University of Hamburg, Germany
    %    email:       glaescher@uke.uni-hamburg.de
    %
    %    $Revision: 1.1 $ $Date: 2004/07/22 09:02:27 $
    
    if isempty(x)
      y = NaN;
      return
    end
    
    if nargin < 2
      dim = find(size(x)~=1, 1 );
      if isempty(dim)
        dim = 1;
      end
    end
    
    
    % Find NaNs in x and nanmean(x)
    nans = isnan(x);
    
    count = size(x,dim) - sum(nans,dim);
    
    
    % Protect against a  all NaNs in one dimension
    i = find(count==0);
    count(i) = 1;
    
    y = nanstd(x,dim)./sqrt(count);
    
    y(i) = i + NaN;
    
    % $Id: nansem.m,v 1.1 2004/07/22 09:02:27 glaescher Exp glaescher $
    
  end

  function Y=create_lut(X)
    
    B = X(end-1:-1:1);
    V = zeros(size(B));
    
    for i = 1:length(B)
      if (i==1)
        V(i) = 1;
      else
        V(i) = prod(B(1:i-1));
      end
    end
    
    Y=zeros(prod(B), length(B));
    
    for i = 1:prod(B)
      
      c = i-1;
      
      for q = length(B):-1:1
        
        a = floor(c/V(q));
        Y(i, length(B)-q+1) = a;
        c = c - a * V(q);
        
      end
    end
    Y = Y + 1;
  end

  function i = isequalcell(c)
    % checks if all cells of a cell array are identical. Works with string and
    % numerical inputs
    
    i = 1;
    c1 = c{1};
    i = 1;
    
    for f=2:length(c)
      
      if ischar(c1)
        if ~strcmp(c1,c{f})
          i = 0;
          break
        end
      elseif isnumeric(c1)
        if c1 ~= c{f}
          i = 0;
          break
        end
      end
      
    end
    
    
  end

end

