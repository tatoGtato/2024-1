function h = axes0(varargin)
% Easily plot the x and y axes through the origin
%
% NOTE: Most of the functionality of axes0 was built directly into the axis 
% object in MATLAB 2015b. See:
% http://www.mathworks.com/help/matlab/creating_plots/display-axis-lines-through-origin.html
%
%% SYNTAX:
%     axes0
%     axes0('PropertyName',propertyvalue,...)
%     h = axes0(...);
%
%% PROPERTIES: (Name-value pairs)
%     'origin' - Position of the origin label
%                'bottomleft'(default)|'topleft'|'topright'|'bottomright'
%                'southwest'|'northwest'|'northeast'|'southeast'|'none'
%     'XTickLabelPos' - Position of the x-axis tick labels with respect to the axis
%                       
%                 'bottom'(default)|'top'
%     'YTickLabelPos' - Position of the y-axis tick labels with respect to the axis                   
%                 'left'(default)|'right'
%     'XLabelPos' - Position of the x-axis label with respect to the axis
%                   'topright (default)'|'bottomright'|'topleft'|'bottomleft'
%                   'firstquad'|'secondquad'|'thirdquad'|'fourthquad'
%                   'left'|'right'
%     'YLabelPos' - Position of the y-axis label with respect to the axis
%                   'topright (default)'|'bottomright'|'topleft'|'bottomleft'
%                   'firstquad'|'secondquad'|'thirdquad'|'fourthquad'
%                   'top'|'bottom'
%      'ticks' - Specify whether to plot the tick markers
%               'on'(default)|'off'|logical
%
%   Legacy properties (still work but have been replaced):
%     'xlabels' - Position of the x-axis labels with respect to the axis
%                 'bottom'(default)|'top'
%     'ylabels' - Position of the y-axis labels with respect to the axis
%                 'left'(default)|'right'

%% DESCRIPTION:
%     h = axes0('PropertyName',propertyvalue,...)
%         Removes the current axes and plots new axes at the origin. Tick
%         labels and markers are preserved by default. Can optionally
%         specify whether to remove tick markers. Can specify the position
%         of x, y, and origin labels. If an output is specified, the
%         function returns handles in the structure h. The following fields
%         are returned:
%           xt = handle of the x tick labels
%           yt = handle of the y tick labels
%           o  = handle of the origin label
%           xl = handle of the x axis label
%           yl = handle of the y axis label
%         All handles are text objects and any text object properties can
%         be modified using set(h.field,'PropertyName',propertyvalue) after
%         axes0 has been called.
%         
%
%% EXAMPLES:
%% Example 1 - Default settings:

% figure;
% x = linspace(-2,2,101);
% plot(x,2*x.^3-3*x+1);
% xlabel('x')
% ylabel('y','Rotation',0)
% 
% axes0
% 
% %% Example 2 - Some Name-Value pairs for better aesthetics: 
% 
% figure;
% plot(x, x.^2-2*x-1);
% xlabel('t')
% ylabel('f(t)','Rotation',0)
% axes0('origin','topright','xticklabelpos','top','yticklabelpos','right',...
%       'xlabelpos','fourthquad','ylabelpos','secondquad');
% 
% %% Example 3 - Modify labels using text object properties:
% 
% figure
% x = linspace(-2,2,101);
% plot(x,2*x.^3-3*x+1);
% h = axes0;
% set(h.xt,'FontSize',15,'color','b')
% set(h.o,'FontSize',15,'color','r','fontweight','bold')
% 
% %% Example 4 - axes0 inherits most properties of the original axis
% 
% figure;
% x = linspace(-1,2,101);
% 
% plot(x,x)
% xlabel('x','fontsize',16)
% ylabel('$y(x)$','Rotation',0,'fontsize',19,'color','b','interpreter','latex')
% set(gca,'fontsize',11,'fontweight','bold','fontname','times','fontangle','italic')
% 
% axes0
%

% Written by Delyle Polet
% dtpolet@ucalgary.ca
% 2015-12-04
% 
% Last update:
% 2016-04-13
%
% The following post on stack exchange inspired this function:
% http://stackoverflow.com/a/2945107/4941405
%
% Copyright (c) 2015, Delyle Polet
% All rights reserved.

% check number of outputs and issue some helpful information
if nargout > 1
    error(['Too many outputs', char(10),...
'If you expected 2 or more outputs, note that handles are now contained in',char(10),...
'a single structure. See documentation for details.'])
end


% Parse inputs

p = inputParser;
addParameter(p, 'ticks', 'on')
validpositions = {'topleft','topright','bottomleft','bottomright',...
                  'northwest','northeast','southwest','southeast','none'};
addParameter(p,'origin', 'bottomleft', @(x) any(validatestring(x,validpositions)))
validpositions = {'topleft','topright','bottomleft','bottomright',...
                  'firstquad','secondquad','thirdquad','fourthquad'};
addParameter(p,'XLabelPos', 'topright', @(x) any(validatestring(x,[validpositions,'left','right'])))
addParameter(p,'YLabelPos', 'topright', @(x) any(validatestring(x,[validpositions,'top','bottom'])))

addParameter(p,'XTickLabelPos','bottom',@(x) any(validatestring(x,{'top','bottom'})))
addParameter(p,'YTickLabelPos','left',@(x) any(validatestring(x,{'left','right'})))

% these parameters kept for legacy
addParameter(p,'xlabels','bottom', @(x) any(validatestring(x,{'top','bottom'})))
addParameter(p,'ylabels','left',@(x) any(validatestring(x,{'left','right'})))

parse(p,varargin{:})

pos.o = p.Results.origin;
pos.x = p.Results.xlabels;
pos.y = p.Results.ylabels;
pos.xl = p.Results.XLabelPos;
pos.yl = p.Results.YLabelPos;

% Overwrite pos.x and pos.y if XTickLabelPos or YTickLabelPos have been
% specified
if ~any(strcmpi(p.UsingDefaults,'XTickLabelPos'))
    pos.x = p.Results.XTickLabelPos;
end
if ~any(strcmpi(p.UsingDefaults,'YTickLabelPos'))
    pos.y = p.Results.YTickLabelPos;
end

% Give the user a warning if using xlabels, ylabels
if ~any(ismember({'xlabels','ylabels'},p.UsingDefaults))
warning(...
['The ''xlabels'' and ''ylabels'' properties will be replaced in an',char(10),...
'upcoming release. Use ''XTickLabelPos'' and ''YTickLabelPos'' respectively'])
end

ticks = p.Results.ticks;
if strcmpi(ticks,'on')
    ticks = true;
elseif strcmpi(ticks,'off')
    ticks = false;
elseif ~islogical(ticks)
    error(['The parameter ''ticks'' must have assignment ''on'', ''off'' or',...
            char(10),'a logical.'])
end

xl = xlim;
yl = ylim;

xlh = get(gca,'XLabel');
ylh = get(gca,'YLabel');
xt = get(gca,'XTick');
yt = get(gca,'YTick');
xtl = get(gca,'XTickLabel');
ytl = get(gca,'YTickLabel');
tl = get(gca,'TickLength'); tl = tl(1);
szf = get(gcf,'Position');
szax = get(gca,'Position');
szax = szax.*szf; % Determine the screen size of the plot area


% figure out which axis is longest- this is important for determining the
% size of the tick markers

[~, I] = max(szax(3:4));


% remove the 0 tick label, otherwise it looks messy
xtl{strcmp(xtl,'0')} = ''; 
ytl{strcmp(ytl,'0')} = ''; 

hold on
plot(xl,[0 0],'k-')
plot([0 0],yl,'k-')
if ticks
    if I == 1
        for i = 1:length(xt)
            patch('xdata',xt([i i]),'ydata',tl*[-1 1]*diff(yl)*szax(3)/szax(4), 'edgecolor','k')
        end
        for i = 1:length(yt)
            patch('xdata',tl*[-1 1]*diff(xl),'ydata',yt([i i]), 'edgecolor','k')
        end
    else
        % must be 2
        for i = 1:length(xt)
            patch('xdata',xt([i i]),'ydata',tl*[-1 1]*diff(yl), 'edgecolor','k')
        end
        for i = 1:length(yt)
            patch('xdata',tl*[-1 1]*diff(xl)*szax(4)/szax(3),'ydata',yt([i i]), 'edgecolor','k')
        end
    end
else
end
xlim(xl)
ylim(yl)


ax = gca;

% Plot tick Labels
h = plotTickLabels(pos,xt,yt,xl,yl,xtl,ytl,tl,I,szax,ax);

% place x and y axis labels
h = placeAxisLabels(xlh,ylh,pos,xl,yl,h,tl);


ax.XLabel.String = '';
ax.YLabel.String = '';
set(gca,'XTick',[], 'XColor',[1 1 1 1])
set(gca,'YTick',[],'Ycolor',[1 1 1 1])
box off

if nargout == 0
    clearvars h
end

end

function h = plotTickLabels(pos,xt,yt,xl,yl,xtl,ytl,tl,I,szax,axh)
    if I == 1
       x_offset = tl+0.005;
       y_offset = (tl+0.005)*szax(3)/szax(4);
    else
        %must be 2
       x_offset = (tl+0.005)*szax(4)/szax(3);
       y_offset = tl+0.005;
    end

    % Plot the x-axis and return a handle for the labels
    if strcmpi(pos.x,'bottom')
        h.xt = text((xt-xl(1))./diff(xl),ones(size(xt)).*(-yl(1)/diff(yl))-y_offset,xtl,...
            'units','normalized','HorizontalAlignment','center','VerticalAlignment','top');
    else
        % must be top
        h.xt = text((xt-xl(1))./diff(xl),ones(size(xt)).*(-yl(1)/diff(yl))+y_offset,xtl,...
            'units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom');      
    end
    textoptions_tick_labels(h.xt,axh);
    
    % Plot the y-axis and return a handle for the labels
    if strcmpi(pos.y,'left')
        h.yt = text(ones(size(yt)).*(-xl(1)/diff(xl))-x_offset,(yt-yl(1))./diff(yl),ytl,...
            'units','normalized','HorizontalAlignment','right','VerticalAlignment','middle');
    else
        % must be right
        h.yt = text(ones(size(yt)).*(-xl(1)/diff(xl))+x_offset,(yt-yl(1))./diff(yl),ytl,...
            'units','normalized','HorizontalAlignment','left','VerticalAlignment','middle');
    end
    textoptions_tick_labels(h.yt,axh);
    
    % Plot the (0,0) marker specially
    if any(strcmpi(pos.o,{'bottomleft','southwest'}))
        h.o = text(-xl(1)/diff(xl)-x_offset,-yl(1)/diff(yl)-y_offset,'0','units','normalized',...
            'VerticalAlignment','top','horizontalalignment','right');
    elseif any(strcmpi(pos.o,{'bottomright','southeast'}))
        h.o = text(-xl(1)/diff(xl)+x_offset,-yl(1)/diff(yl)-y_offset,'0','units','normalized',...
            'VerticalAlignment','top');
    elseif any(strcmpi(pos.o,{'topleft','northwest'}))
        h.o = text(-xl(1)/diff(xl)-x_offset,-yl(1)/diff(yl)+y_offset,'0','units','normalized',...
            'VerticalAlignment','bottom','horizontalalignment','right');
    elseif any(strcmpi(pos.o,{'topright','northeast'}))
        h.o = text(-xl(1)/diff(xl)+x_offset,-yl(1)/diff(yl)+y_offset,'0','units','normalized',...
            'VerticalAlignment','bottom');
    else
        % If 'none' specified, an empty label is plotted in the southwest
        % corner
        h.o = text(-xl(1)/diff(xl)-x_offset,-yl(1)/diff(yl)-y_offset,'','units','normalized',...
            'VerticalAlignment','top','horizontalalignment','right');
    end
    textoptions_tick_labels(h.o,axh);
    
end

function h = placeAxisLabels(xlh,ylh,pos,xl,yl,h,tl)
    % set xlabel initially to topright (the default)
    ml = D2N(0,yl); % y position of the xaxis
    hxl = text(1.01,ml,xlh.String,'units','normalized','verticalalignment','bottom');
    textoptions_axis_labels(hxl,xlh);
    
    xt_extent = h.xt.Extent;
    if strcmpi(pos.xl,'bottomright')
        set(hxl,'verticalalignment','top');
    elseif strcmpi(pos.xl,'topleft')
        set(hxl,'position',[-0.01,ml],'horizontalalignment','right','verticalalignment','bottom');
    elseif strcmpi(pos.xl,'bottomleft')
        set(hxl,'position',[-0.01,ml],'horizontalalignment','right','verticalalignment','top');
    elseif strcmpi(pos.xl,'right')
        set(hxl,'verticalalignment','middle');
    elseif strcmpi(pos.xl,'left')
        set(hxl,'position',[-0.01,ml],'horizontalalignment','right','verticalalignment','middle');
    elseif strcmpi(pos.xl,'firstquad')
        set(hxl,'Position',[D2N(xl(2)/2,xl),max(sum(xt_extent([2,4]))+0.005,tl/2+ml+0.02)],'horizontalAlignment','center','verticalalignment','bottom')
    elseif strcmpi(pos.xl,'secondquad')
        set(hxl,'Position',[D2N(xl(1)/2,xl),max(sum(xt_extent([2,4]))+0.005,tl/2+ml+0.02)],'horizontalAlignment','center','verticalalignment','bottom')
    elseif strcmpi(pos.xl,'thirdquad')
        set(hxl,'Position',[D2N(xl(1)/2,xl),min(xt_extent(2)-0.005,-tl/2+ml-0.02)],'horizontalAlignment','center','verticalalignment','top')
    elseif strcmpi(pos.xl,'fourthquad')
        set(hxl,'Position',[D2N(xl(2)/2,xl),min(xt_extent(2)-0.005,-tl/2+ml-0.02)],'horizontalAlignment','center','verticalalignment','top')
    end
    
    % set xlabel initially to topright (the default)
    ml = D2N(0,xl); % x position of the y axis
    hyl = text(ml,1.01,ylh.String,'units','normalized','verticalalignment','bottom');
    textoptions_axis_labels(hyl,ylh);
    mxyt = max_extent(h.yt);
    if strcmpi(pos.yl,'bottomright')
        set(hyl,'position',[ml,-0.01],'verticalalignment','top');
    elseif strcmpi(pos.yl,'topleft')
        set(hyl,'horizontalalignment','right','verticalalignment','bottom');
    elseif strcmpi(pos.yl,'bottomleft')
        set(hyl,'position',[ml,-0.01],'horizontalalignment','right','verticalalignment','top');
    elseif strcmpi(pos.yl,'top')
        set(hyl,'verticalalignment','bottom','horizontalalignment','center');
    elseif strcmpi(pos.yl,'bottom')
        set(hyl,'position',[ml,-0.01],'horizontalalignment','center','verticalalignment','top');
    elseif strcmpi(pos.yl,'firstquad')
        set(hyl,'Position',[max(mxyt.right+0.005,tl/2+ml+0.02),D2N(yl(2)/2,yl)],'horizontalAlignment','left','verticalalignment','middle')
    elseif strcmpi(pos.yl,'secondquad')
        set(hyl,'Position',[min(mxyt.left-0.005,-tl/2+ml-0.02),D2N(yl(2)/2,yl)],'horizontalAlignment','right','verticalalignment','middle')
    elseif strcmpi(pos.yl,'thirdquad')
        set(hyl,'Position',[min(mxyt.left-0.005,-tl/2+ml-0.02),D2N(yl(1)/2,yl)],'horizontalAlignment','right','verticalalignment','middle')
    elseif strcmpi(pos.yl,'fourthquad')
        set(hyl,'Position',[max(mxyt.right+0.005,tl/2+ml+0.02),D2N(yl(1)/2,yl)],'horizontalAlignment','left','verticalalignment','middle')
    end
    h.xl = hxl;
    h.yl = hyl;
end

function textoptions_axis_labels(th,ih)
% inherit most of the text options from the initial handle ih
options = {'FontSize','FontWeight','FontName','Color','BackgroundColor',...
    'EdgeColor','FontAngle','FontSmoothing','FontUnits','Interpreter',...
    'LineStyle','LineWidth','Rotation'};

for i = 1:length(options)
    th.(options{i}) = ih.(options{i});
end
end

function textoptions_tick_labels(th,ih)
% inherit most of the text options from the initial handle ih
options = {'FontSize','FontWeight','FontName',...
    'FontAngle','FontSmoothing','FontUnits'};
for j = 1:size(th,1)
    for i = 1:length(options)
        th(j).(options{i}) = ih.(options{i});
    end
end
end

function n = D2N(d,lim)
% converts from data units to normalized units
n = (d-lim(1))/diff(lim);
end

function mx = max_extent(text_array)
    M = NaN(size(text_array,1),4);
    for i = 1:size(text_array,1)
         M(i,:) = text_array(i).Extent;
    end
    mx.left = min(M(:,1));
    mx.right = max(sum(M(:,[1,3]),2));
    mx.top = max(sum(M(:,[2,4]),2));
    mx.bottom = max(M(:,2));
end