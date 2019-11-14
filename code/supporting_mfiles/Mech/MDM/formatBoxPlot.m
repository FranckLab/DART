function [h, a] = formatBoxPlot(wh)
if nargin < 1, wh = [1 1]; end


% wh = [1 1];
h = gcf;

set(gca,'color','none')
set(h,'color','none');
set(h, 'PaperUnits','inches')
set(h, 'PaperSize',wh)
set(h, 'Paperposition', [3 3 wh]);
set(h, 'Units','inches','Position', [3, 3, wh]);



a = gca;
set(a,'units','points')
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(a,'fontsize',8,'fontname','arial');
set(a,'Linewidth',1,'color','none')
set(a, 'LooseInset', get(a,'TightInset'));
set(a,'xtick', []);




set(findobj(gcf,'Type','text'),'FontSize',8,'fontname','arial')
box off


set(gca,'color','none')
set(h,'color','none');
% set(h, 'InvertHardCopy', 'off')
set(h, 'PaperUnits','inches')
set(h, 'PaperSize',wh)
set(h, 'Paperposition', [3 3 wh]);
set(h, 'Units','inches','Position', [3, 3, wh]);
        
        
        
% set(gcf,'papersize',[3.25,2],'position',[6 4 3.25 2]);