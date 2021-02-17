


f1=figure('color','w','position',[100 100 1000 800]);
ax11 = axes('position',[.1 .70 .74 .28]); hold on; box on; grid on;
ax12 = axes('position',[.1 .39 .74 .28]); hold on; box on; grid on;
ax13 = axes('position',[.1 .08 .74 .28]); hold on; box on; grid on;
set([ax11 ax12 ax13],'ticklabelinterpreter','latex','fontsize',18,'xlim',...
    [1 max(GPeaks(:,9))])
set([ax11,ax12],'xticklabel',[])

% Positions. 
p11 = plot(ax11,GPeaks(:,9),GPeaks(:,2),'.','markersize',12,'color',cm(1,:)); hold on
p12 = plot(ax11,GPeaks(:,9),GPeaks(:,5),'o','color',cm(1,:));
p13 = plot(ax11,GPeaks(:,9),GPeaks(:,8),'x','markersize',10,'color',cm(1,:));
ylim= get(ax11,'ylim');
yticks = ylim(1):.001:ylim(2)+.001;
for mm = 1:length(yticks)
    plot(ax11,[1 max(GPeaks(:,9))],[1 1]*yticks(mm),'linewidth',.5,'color',[.5 .5 .5 .25])
end
leg = legend(ax11,[p11 p12 p13],'FindPeaks','Gfit','MainPeak');
set(leg,'interpreter','latex','box','on','position',[.87 .85 .1 .1]);

% Ampltidues. 
p21 = plot(ax12,GPeaks(:,9),GPeaks(:,1),'.','markersize',12,'color',cm(1,:)); hold on
p22 = plot(ax12,GPeaks(:,9),GPeaks(:,4),'o','color',cm(1,:));
p23 = plot(ax12,GPeaks(:,9),GPeaks(:,7),'x','markersize',10,'color',cm(1,:));
% leg = legend(ax2,[p11 p12 p13],'MaxP','Gfit','MaxMean');
% set(leg,'interpreter','latex','box','off');
p24 = plot(ax12,[0 max(GPeaks(:,9))],[1 1]*482,'-','color',cm(2,:)); % saturation
text(ax12,1.01,.8,'saturation','interpreter','latex','units','normalized',...
    'color',cm(2,:),'fontsize',14)
text(ax12,1.01,1.4,['$\Delta t=$',num2str( abs(tt(1,1)-tt(1,2)) ),' $\mu$s',...
    ''],...
    'interpreter','latex','units','normalized',...
    'color',[1 1 1]*.5,'fontsize',14)
% text(ax2,1.01,1.3,[num2str( round(1e6*2*Da*e*( (0.001/1e6)/d *sqrt(E0*s0+E1*s1) ).^2,2) ),'e-6 u'
% ],...
%     'interpreter','latex','units','normalized',...
%     'color',cm(5,:),'fontsize',14)

% Widths. 
p31 = plot(ax13,GPeaks(:,9),GPeaks(:,3),'.','markersize',12,'color',cm(1,:)); hold on
p32 = plot(ax13,GPeaks(:,9),GPeaks(:,6),'o','color',cm(1,:));
% leg = legend(ax3,[p11 p12],'MaxP','Gfit');
% set(leg,'interpreter','latex','box','off');

xlabel(ax13,'shot no','interpreter','latex')
ylabel(ax11,'position (m/z)','interpreter','latex')
ylabel(ax12,'amplitude (V)','interpreter','latex')
ylabel(ax13,'widths (m/z)','interpreter','latex')



h   = gcf;
fig_name = strrep([wd,'Positions_Amplitudes_Win',num2str(round(min(xwin),2)),'-',...
        num2str(round(max(xwin),2)),'us'],'.',',');
set(h,'unit','inches');
pos = get(h,'position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
print(strcat(fig_name),'-dpng')


% close(f1)