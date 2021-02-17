
clear variables
%__________________________________________________________________________
%__________________________________________________________________________
%##########################################################################
% Type in the path, where the results are stored. 
wd      = 'G:\Gabriele\Raw_Data\Results\2019-10-08\';
%##########################################################################


%__________________________________________________________________________
% Each textfile refers to a specific window, where peaks were evaluated.
% The main peak (=maximum peak) is saved into the textfile for a whole
% measurement. The maximum peak position can deviate drastically from each 
% each other if a specific shot did not contain the peak defined as main 
% peak in the rest of the image series. 
% Each matfile contains all peaks found in a single shot.  
txtfiles    = dir([wd,'*.txt']); txtfiles={txtfiles.name}';sigfiles=txtfiles;
matfiles    = dir([wd,'*.mat']); matfiles={matfiles.name}';
TF1         = contains(txtfiles,'MaxPeak'); 
TF2         = contains(txtfiles,'SigPeaks');
sigfiles    = sigfiles(TF2==1);txtfiles=txtfiles(TF1==1);

% Read out time window of each file. 
filesall    = {txtfiles;sigfiles;matfiles}; numall = filesall;
for jj      = 1:length(filesall)
    filenames   = filesall{jj};
    num         = zeros(length(filenames),2); num_sort = 1;
    for ii      = 1:length(filenames)
        num_ii  = regexp(strrep(filenames{ii},',','.'),'(+|-)?\d+(\.\d+)?(E(+|-)?\d+)?','Match');%\d*
        num(ii,:)= abs( [str2double(num_ii{1}) str2double(num_ii{end})] );
    end
    num0        = num(:,1)+abs(num(:,1)-num(:,2)); % t of the time window
    numall{jj}  = [num0 num];
end
%__________________________________________________________________________


 


%% Consider the main peak in a specific window (of about 1 proton mass width)
% Goal is to plot the single shot peak position and amplitude to see the
% shot-to-shot fluctuation.
%##########################################################################
% Pick 3 windows (=files) you want to compare. Window at the raw us scale. 
w1      = .5; % [us]
w2      = 3.4; % [us]
w3      = 6.; % [us]
wall    = [w1 w2 w3];
%##########################################################################

%__________________________________________________________________________
xsz = .35; xsh = .44;
figure('color','w','position',[100 100 1000 800]); cm = colormap('lines');
ax1 = axes('position',[.1 .70 xsz .28]); hold on; box on; grid on;
ax2 = axes('position',[.1 .39 xsz .28]); hold on; box on; grid on;
ax3 = axes('position',[.1 .08 xsz .28]); hold on; box on; grid on;
ax4 = axes('position',[.1+xsh .70 xsz .28]); hold on; box on; grid on;
ax5 = axes('position',[.1+xsh .39 xsz .28]); hold on; box on; grid on;
ax6 = axes('position',[.1+xsh .08 xsz .28]); hold on; box on; grid on;
set([ax1 ax2 ax3 ax4 ax5 ax6],'ticklabelinterpreter','latex','fontsize',18)
set([ax1,ax2,ax4,ax5],'xticklabel',[])
axall = [ax1 ax2 ax3 ax4 ax5 ax6];
set([ax4 ax5 ax6],'yaxislocation','right')
xlabel(ax3,'shot no','interpreter','latex')
ylabel(ax2,'position (m/z)','interpreter','latex')
xlabel(ax6,'amplitude','interpreter','latex')
ylabel(ax5,'position (m/z)','interpreter','latex')
%__________________________________________________________________________

for ii          = 1:length(wall)
    num         = numall{1}; num0 = num(:,1); num = num(:,2:3);
    [~,posw]    = min(abs(num0-wall(ii)));
    TF          = contains(txtfiles,strrep(num2str(num(posw,1)),'.',','));
    txtfiles1   = txtfiles(TF);
    GPeaks      = load([wd,txtfiles1{1}]);
    % [MaxPksAmp MaxPksPos MaxPksWdt GPksAmp GPksPos GPksWdt MaxPmeanAmp MaxPPos ShotNo] 
    
    % Positions. 
    p11 = plot(axall(ii),GPeaks(:,9),GPeaks(:,2),'.','markersize',12,'color',cm(1,:)); hold on
    p12 = plot(axall(ii),GPeaks(:,9),GPeaks(:,5),'o','color',cm(1,:));
    p13 = plot(axall(ii),GPeaks(:,9),GPeaks(:,8),'x','markersize',10,'color',cm(1,:));
    % Draw the pixel distance. 
    ylim= get(axall(ii),'ylim');
    yticks = ylim(1):.001:ylim(2)+.001;
    for mm = 1:length(yticks)
        plot(axall(ii),[1 max(GPeaks(:,9))],[1 1]*yticks(mm),'linewidth',.5,'color',[.5 .5 .5 .25])
    end
    % Amplitude.
    yyaxis(axall(ii),'right')
    plot(axall(ii),GPeaks(:,9),GPeaks(:,7),'.','markersize',12,'color',cm(5,:)); hold on
    plot(axall(ii),[1 max(GPeaks(:,9))],[1 1]*481,'color',cm(2,:))
    set(axall(ii),'ycolor',cm(5,:))
    if ii==2
        ylabel(axall(ii),'amplitude (m/z)','interpreter','latex')
    end
    yyaxis(axall(ii),'left')
    
    %**********************************************************************
    LThresh = 100;
    UThresh = 400;
    p21 = plot(axall(ii+3),GPeaks(:,7),GPeaks(:,2),'.','markersize',12,'color',cm(1,:)); hold on
    p22 = plot(axall(ii+3),GPeaks(:,7),GPeaks(:,5),'o','color',cm(1,:));
    p23 = plot(axall(ii+3),GPeaks(:,7),GPeaks(:,8),'x','markersize',10,'color',cm(1,:));
    % Mean value and deviation.
    vari= 2;
    plot(axall(ii+3),[1 1]*LThresh,[ylim(1) ylim(2)],'--','color',cm(2,:))
    plot(axall(ii+3),[1 1]*UThresh,[ylim(1) ylim(2)],'--','color',cm(2,:))
    TF  = and(GPeaks(:,7)>LThresh,GPeaks(:,7)<UThresh);
    xpts= GPeaks(TF,9);
    ypts= GPeaks(TF,vari);
    ymed= median(ypts,'omitnan');
    ymad= mad(ypts);
    plot(axall(ii+3),[min(GPeaks(:,7)) 500],[1 1]*ymed,'--','color',[.5 .5 .5 .5],'linewidth',2)
    disp(ymad)
    patch(axall(ii+3),[min(GPeaks(:,7)) 500 500,...
        min(GPeaks(:,7)) min(GPeaks(:,7))],...
        ymed+[-ymad -ymad ymad ymad -ymad],[.5 .5 .5],'edgecolor','none',...
        'facealpha',.5)
    text(axall(ii+3),.01,.05,['shot-to-shot-fluct: ',num2str(round(ymed,3)),...
        '$\pm$',num2str(round(ymad,3)),' $\mu$s'],...
        'units','normalized','interpreter','latex','fontsize',14)
    %**********************************************************************
end

leg = legend(ax4,[p11 p12 p13],'FindPeaks','Gfit','MainPeak');
set(leg,'interpreter','latex','box','on','location','best');%'position',[.87 .85 .1 .1]);

set(axall(1:3),'xlim', [1 max(GPeaks(:,9))])
set(axall(4:6),'xlim', [1 500])





h           = gcf;
fig_name    = strrep([wd,'MaxSignalPeakPlot_',num2str(wall),'us'],'.','');
set(h,'unit','inches');
pos         = get(h,'position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
print(strcat(fig_name),'-dpng')


%% Consider all peaks over threshold in a specific window (of about 1 proton mass width)
%##########################################################################
% Pick 1 window (=files). Window at the raw us scale. 
w1      = wall(1); % [us]
% Pick 3 intervals where to find a mean value.
thresh  = [ .64 .66;...
            .66 .68;...
            .69 .71];
AmplThresh = 470;
%##########################################################################

%__________________________________________________________________________
xsz = .85; xsh = .44;
figure('color','w','position',[100 100 1000 800]); cm = colormap('lines');
ax1 = axes('position',[.1 .70 xsz .28]); hold on; box on; grid on;
ax2 = axes('position',[.1 .39 xsz .28]); hold on; box on; grid on;
ax3 = axes('position',[.1 .08 xsz .28]); hold on; box on; grid on;
% ax4 = axes('position',[.1+xsh .70 xsz .28]); hold on; box on; grid on;
% ax5 = axes('position',[.1+xsh .39 xsz .28]); hold on; box on; grid on;
% ax6 = axes('position',[.1+xsh .08 xsz .28]); hold on; box on; grid on;
set([ax1 ax2 ax3],'ticklabelinterpreter','latex','fontsize',18)% ax4 ax5 ax6
set([ax1,ax2],'xticklabel',[])%ax4,ax5
axall = [ax1 ax2 ax3];% ax4 ax5 ax6
% set([ax4 ax5 ax6],'yaxislocation','right')
xlabel(ax3,'amplitude','interpreter','latex')
ylabel(ax2,'position (m/z)','interpreter','latex')
%__________________________________________________________________________

num         = numall{2}; num0 = num(:,1); num = num(:,2:3);
[~,posw]    = min(abs(num0-w1));
TF          = contains(sigfiles,strrep(num2str(num(posw,1)),'.',','));
sigfiles1   = sigfiles(TF); 
SPeaks      = load([wd,sigfiles1{1}]); % row 1 is shot number, first half is position, second half amplitude
SPeaks(SPeaks==0)=nan;
PS2         = (size(SPeaks,2)-2)/2; % number of peaks above threshold

%__________________________________________________________________________
for ii = 1:3
    vari    = PS2+3:size(SPeaks,2); % amplitude values
    % Positions.
    for kk  = 2:PS2+1
        p11 = plot(axall(ii),SPeaks(:,kk+PS2+1),SPeaks(:,kk),'.','markersize',12,'color',cm(kk,:)); hold on
    end
    % Draw the pixel distance.
    ylim    = get(axall(ii),'ylim');
    yticks  = ylim(1):.001:ylim(2)+.001;
    for mm  = 1:length(yticks)
        plot(axall(ii),[1 max(SPeaks(:,1))],[1 1]*yticks(mm),'linewidth',.5,'color',[.5 .5 .5 .25])
    end
    
    %**********************************************************************
    LThresh = thresh(ii,1); % position threshold
    UThresh = thresh(ii,2); % position threshold
    % Mean value and deviation.
    p1 = patch(axall(ii),[min(min(SPeaks(:,vari))) 500 ...
        500 min(min(SPeaks(:,vari))) min(min(SPeaks(:,vari)))],...
        [[1 1]*LThresh [1 1]*UThresh LThresh],'k','edgecolor','none',...
        'facealpha',.2);
    xlim = max(max(SPeaks(:,vari)));
    patch(axall(ii),[[1 1]*AmplThresh xlim xlim AmplThresh],...
        [ylim(1) ylim(2) ylim(2) ylim(1) ylim(1)],cm(2,:),'edgecolor','none',...
        'facealpha',.5)
    TF      = SPeaks(:,vari)<=AmplThresh;  % count only non-saturated peaks
    xpts    = SPeaks(:,PS2+3:end);    % amplitudes
    ypts    = SPeaks(:,2:PS2+1);    % positions
    xpts(TF==0)=nan;ypts(TF==0)=nan;
    TF      = and(ypts>LThresh,ypts<UThresh);
    xpts(TF==0)=nan;ypts(TF==0)=nan;
    pin     = plot(axall(ii),xpts,ypts,'o','color','k');
    ymed    = median(ypts(:),'omitnan');
    ymad    = mad(ypts(:));
    plot(axall(ii),[1 max(max(SPeaks(:,vari)))],[1 1]*ymed,'-','color',[.5 .5 .5 .5],'linewidth',2)
    patch(axall(ii),[1 max(max(SPeaks(:,vari))) max(max(SPeaks(:,vari))),...
        1 1],...
        ymed+[-ymad -ymad ymad ymad -ymad],[.5 .5 .5],'edgecolor','none',...
        'facealpha',.5)
    text(axall(ii),.01,.05,['shot-to-shot-fluct: ',num2str(round(ymed,3)),...
        '$\pm$',num2str(round(ymad,3)),' $\mu$s'],...
        'units','normalized','interpreter','latex','fontsize',14)
    %**********************************************************************
    
end
%__________________________________________________________________________
% leg = legend(ax2,[p11 p12 p13],'SigAbvThresh','Gfit','MainPeak');
% set(leg,'interpreter','latex','box','on','location','best');%'position',[.87 .85 .1 .1]);
set(axall(1:3),'xlim', [min(min(SPeaks(:,vari))) max(max(SPeaks(:,vari)))])
%__________________________________________________________________________
% Save result as a figure.
h           = gcf;
fig_name    = strrep([wd,'AllSignalAboveThreshPlot_',num2str(wall),'us'],'.','');
set(h,'unit','inches');
pos         = get(h,'position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
print(strcat(fig_name),'-dpng')
%__________________________________________________________________________
   
    
    
    
