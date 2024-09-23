function [ax,h]=master_plot_title(text,fontsize,topleft)
%
%Centers a title over a group of subplots.
%Returns a handle to the title and the handle to an axis.
% [ax,h]=subtitle(text)
%           returns handles to both the axis and the title.
% ax=subtitle(text)
%           returns a handle to the axis only.
if (nargin < 3), topleft = 0; end
if (nargin < 2), fontsize = 16; end
if (~topleft)
    ax=axes('Units','Normal','Position',[.075 .075 .85 .88],'Visible','off');
else
%    ax=axes('Units','Normal','Position',[.075 .075 .15 .85],'Visible','off');
    ax=axes('Units','Normal','Position',[.06 .06 .10 .90],'Visible','off');
end
set(get(ax,'Title'),'Visible','on')
title(text,'fontweight','bold','fontsize',fontsize,'Interpreter','none');

if (nargout < 2)
    return
end
h=get(ax,'Title');