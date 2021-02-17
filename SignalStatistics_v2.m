% Fit sum of Gaussians to a signal (which is not a sum of Gaussians) but to
% understand the beaviour of the data.
clear variables
addpath('C:\Users\tauscher\AppData\Local\MaxPlanck\GunSimulation\Matlab\Sonstige\Andrey')
addpath('G:\Gabriele\massspectrumanalysis Backup 2019-08-01')
addpath('G:\Gabriele\massspectrumanalysis Backup 2019-08-01\plotting')
addpath('G:\Gabriele\massspectrumanalysis Backup 2019-08-01\subroutines')

wd          = 'G:\Gabriele\Raw_Data\Results\2019-10-14\';
wd_dat      = 'G:\Gabriele\Raw_Data\';

filename    = [wd_dat,'2019-10-14\02.signal.div'];%'Measurement_08_10_13\08.signal.div';
logfilename = [wd_dat,'2019-10-14\02.log'];% 'Measurement_08_10_13\08.log';
% filename    = [wd_dat,'Measurement_08_10_13\08.signal.div'];%'Measurement_08_10_13\08.signal.div';
% logfilename = [wd_dat,'Measurement_08_10_13\08.log'];%2019-10-08\03.log 'Measurement_08_10_13\08.log';
% zB 132 (oder 10.000) Schuss in einer Messung. 
vari        = 1:1500;
[tt,yy]     = FCT_read_AndreysData(filename,logfilename,vari); % all shots
 

%% Umrechnung TOF zu Masse ************************************************
Da      = 6.0221412927*10^23*1000;
e       = 1.602176462*10^(-19); % [C] electric charge

% Andrey Kalibrierung, eigentlich individual for single shots
d       = 1.050;
s0      = 0.004;
s1      = 0.047;
V0      = 16000;
V1      = 14500;

%__________________________________________________________________________
E0      = (V0-V1)/s0;
E1      = V1/s1;
% M_conv  = 2*Da*e*( tt/d *sqrt(E0*s0+E1*s1) ).^2; % tt ist der x-Messwert
%__________________________________________________________________________
M_max   = 2*Da*e*( (40/1e6)/d *sqrt(E0*s0+E1*s1) ).^2; % 1.007 u
M_x     = 1.007276466621*linspace(0,round(M_max),round(M_max));
M_proton= 2*Da*e*( (0.651/1e6)/d *sqrt(E0*s0+E1*s1) ).^2; % 1.007 u
dxproton= ( sqrt(1.007/(2*Da*e)) *d*1e6 )  *1/sqrt(E0*s0+E1*s1) ; % 0.5997
x_M     = ( sqrt(M_x/(2*Da*e)) *d*1e6 )  *1/sqrt(E0*s0+E1*s1) ; % waere die korrekte Massen x-Achse zu obiger Kalibrierung
x_tt    = 2*Da*e*( (tt(1,:)/1e6)/d *sqrt(E0*s0+E1*s1) ).^2;
x_test  = ( sqrt((1.*M_x+0*M_x(200))/(2*Da*e)) *(d*1.3)*1e6 )  *1/sqrt(E0*s0+E1*s1) ; % waere die korrekte Massen x-Achse zu obiger Kalibrierung
figure,plot(x_test), hold on,plot(x_M)
%__________________________________________________________________________
%__________________________________________________________________________

% Fenster um ganzzahlige Massen +- 1/2 mp.
xwin    = 1:round(max(x_tt));
wInd    = nan*ones(length(xwin),2);
for ii  = 1:length(xwin)
    [~,posW1]   = min(abs(x_tt-(xwin(ii)-.5)));
    [~,posW2]   = min(abs(x_tt-(xwin(ii)+.5)));
   wInd(ii,:)   = [(posW1) (posW2)]; % [x_tt(posW1) x_tt(posW2)]; 
end
%__________________________________________________________________________
%__________________________________________________________________________


%**************************************************************************
%**************************************************************************
% Find peaks in the windowed data and plot corresponding Gaussian.
no_pks      = 3; % ,'Npeaks',no_pks for how many peaks you want to look
min_dist    = dxproton/1000; % minimum peak distance
min_height  = 15; % minimum peak height (mind. 3*noise Faustregel)
min_prom    = 10; % minimum peak prominence (i.e. noise)
thresh      = 15; % unter der threshold ist noise (per definition here) 'Threshold',thresh,
SigThresh   = 200;% alle über der threshold werden ins file geschrieben
        
%**************************************************************************
case_fitequ = 'Gauss'; % ['Gauss','asymmG','sumfit']
wcutR       = 3.; % cutoff of the data to the right of the peak position
wcutL       = 3; %
%**************************************************************************
%**************************************************************************


%%
%__________________________________________________________________________
fCD     = figure('position',[1800 100 1400 800],'color','w'); cm = colormap('lines');
xank    = .1; yank = .73; xsz = .4; ysz = .24; ysh = .08; xsh = .1;
ax1     = axes('position',[xank yank-0*(ysz+ysh) xsz ysz]); hold on; box on
ax2     = axes('position',[xank yank-1*(ysz+ysh) xsz ysz]); hold on; box on
ax3     = axes('position',[xank yank-2*(ysz+ysh) xsz ysz]); hold on; box on
ax4     = axes('position',[xank+1*(xsh+xsz) yank-0*(ysz+ysh) xsz*.8 ysz]); hold on; box on
ax5     = axes('position',[xank+1*(xsh+xsz) yank-1*(ysz+ysh) xsz*.8 ysz]); hold on; box on
ax6     = axes('position',[xank+1*(xsh+xsz) yank-2*(ysz+ysh) xsz*.8 ysz]); hold on; box on
set([ax1 ax2 ax3 ax4 ax5 ax6],'ticklabelinterpreter','latex','fontsize',16)
xlabel(ax3,'$x$ (mass/charge)','interpreter','latex')
ylabel(ax2,'signal amplitude (V)','interpreter','latex')
xlabel(ax4,'noise/residuum (V)','interpreter','latex')
ylabel(ax4,'probability','interpreter','latex')
set(ax4,'yaxislocation','right')

xlabel(ax6,'peak position $x$ (mass/charge)','interpreter','latex')
ylabel(ax5,'peak ampl (V)','interpreter','latex')
ylabel(ax6,'peak width ($x$)','interpreter','latex')

%**************************************************************************
% data output.
%**************************************************************************

for kk          = 100:200%:length(wInd) % move the window from left to right
    get(fCD)
    GPeaks      = nan*ones(1,9);    % [Mampl Mpos Mwidth Gampl Gpos Gwidth Maxpos meanAmpl shotno]
    SPeaksPos   = nan*ones(1,50);   
    SPeaksAmp   = SPeaksPos;
    for jj      = 1:132
        %______________________________________________________________________
        % Load some data with peaks.
        xdat    = tt(jj,:);%2*Da*e*( tt(jj,:)/d *sqrt(E0*s0+E1*s1) ).^2;%tt(jj,:);
        ydat    = yy(jj,:);
        %______________________________________________________________________
        % Sättigung an oder aus?
        ydat0   = ydat; xdat0 = xdat;
        ydat(ydat>=482)=NaN; % gesaettigt
        %______________________________________________________________________
        % Define noise by threshold value. Find noise peaks to define the
        % threshold.
        ydat1   = ydat0;
        ydat1(ydat1>thresh)     = nan;
        baselineEst             = median(ydat1,'omitnan');
%         [pks0,locs0,w0,prom0]   = findpeaks(ydat0,xdat0);
        %______________________________________________________________________
        %______________________________________________________________________
        % Define window, where you want to look at the data.
        xstart  = wInd(kk,1); % start x-value of the considered window (left)
        xstop   = wInd(kk,2); % end x-value of the considered window (right)
        xwin    = xdat(xstart:xstop);
        ywin    = ydat(xstart:xstop);
        %______________________________________________________________________
        % Max peak finder and starting conditions for Gaussian fits. 
        [pks,locs,w,~]          = findpeaks(ywin,xwin,'MinPeakDist',min_dist,...
            'MinPeakHeight',min_height,'MinPeakProminence',min_prom); % ,'MinPeakDistance',min_dist
        if isempty(pks)==1,continue, end
        % For the plot only. Control display. 
        Gauss   = nan*ones(length(locs),length(xwin)); % preallocated dummy
        xx      = linspace(xwin(1),xwin(end),1000);
        for n   = 1:length(locs)
            Gauss(n,:)          =  pks(n)*exp(-((xwin - locs(n))/w(n)).^2);
        end
        %______________________________________________________________________
        switch case_fitequ
            case 'Gauss'
                fiteq           = fittype( 'a1/sqrt(2*pi)/c1*exp(- (x-b1).^2./2./c1^2 )' ); % Gauss fitequation
                Gfits           = cell(length(locs),1);
            case 'asymmG' % hier eigentlich zu wenig Datenpunkte
                fiteq           = fittype( '(a1+d1.*(x-b1)/c1)./( (x-b1).^2/c1.^2+1 )' ); % asymm Gauss fitequation
                Gfits           = cell(length(locs),1);
            case 'sumfit' % asymm Gauss (otherwise change in the beginning)
                %__________________________________________________________________________
                I_lambda_oo     = []; % sum fit equation (je nachdem wie viele peaks es gab) 'F0+F1.*x' Untergrund
                for eq_oo       = 1:length(pks)
                    Ilambda_eq  = ['(a',num2str(eq_oo+9),'+f',num2str(eq_oo+9),'.*(x-b',num2str(eq_oo+9),')/c',num2str(eq_oo+9),')./((x-b',num2str(eq_oo+9),').^2/c',num2str(eq_oo+9),'.^2+1)']; %+(I_l2+(A_l2.*(x-b2)/c2))./((x-b2).^2/c2.^2+1)
                    I_lambda_oo = [I_lambda_oo,'+',Ilambda_eq];
                end
                %__________________________________________________________________________
                fiteq           = I_lambda_oo;
                Gfits           = cell(1,1);
        end
        if isempty(fiteq)==1, continue, end
        FO                      = fitoptions(fiteq); % fitoptions
        %__________________________________________________________________________
        %**************************************************************************
        % Fit a Gaussian (better than the findpeaks function does). If the fit
        % looks bad, tighten the upper and lower bonds in the fitoptions FO. Start
        % conditions are taken from findpeaks.
        Gampl       = nan*ones(length(locs),1);
        Gpos        = nan*ones(length(locs),1);
        Gwds        = nan*ones(length(locs),1);
        GaussfitPlot= zeros(length(xwin),1);
        xyPl        = cell(1,1);
        if strcmp(case_fitequ,'sumfit')
            FO.StartPoint       = [   pks(:)' locs(:)'       2*w(:)'    zeros(length(pks),1)']; % start point [a1 b1 c1], where a1 is the peak amplitude, b1 is the x-position and c1 is the sigma of the Gaussian
            FO.Upper            = [10*pks(:)' locs(:)'+w(:)'   w(:)'*4  ones(length(pks),1)' ]; % upper limits of the coefficients
            FO.Lower            = [ 0*pks(:)' locs(:)'-w(:)'   w(:)'/4 -ones(length(pks),1)' ]; % lower limits of the coefficients
            xfit                = xwin(:);
            yfit                = ywin(:); 
            yfit(isnan(xfit))   = [];xfit(isnan(xfit)) = [];
            xfit(isnan(yfit))   = [];yfit(isnan(yfit)) = [];
            Gaussfit            = fit(xfit(:),yfit(:),fiteq,FO);
            GaussfitPlot        = feval(Gaussfit,xfit);
            Gfits{1}            = Gaussfit;
            xyPl                = [xfit,yfit]; xyPl={xyPl};
            %______________________________________________________________________
            % save position and amplitude of the Gaussianfit:
            Gcoeff              =   coeffvalues(Gaussfit);
            Gampl               = Gcoeff(1:length(locs))';
            Gpos                = Gcoeff(length(locs)+1:length(locs)*2)';
            Gwds                = Gcoeff(length(locs)*2+1:length(locs)*3)';
            %______________________________________________________________________
        else
            for n = 1:length(locs) % each peak is fitted individually
                % For the fit, only consider the peak of this loop n. Or fit the sum
                % instead, is that more accurate?
                [~,xstart_n]    = min(abs( xwin- ( locs(n)-wcutL*w(n) ) ));
                [~,xstop_n]     = min(abs( xwin- ( locs(n)+wcutR*w(n) ) ));
                xfit            = xwin(xstart_n:xstop_n);
                yfit            = ywin(xstart_n:xstop_n);
                yfit(isnan(xfit))= [];xfit(isnan(xfit)) = [];
                xfit(isnan(yfit))= [];yfit(isnan(yfit)) = [];
                switch case_fitequ
                    case 'Gauss'
                        FO.StartPoint   = [  2*pi*pks(n) locs(n)      2*w(n)  ]; % start point [a1 b1 c1], where a1 is the peak amplitude, b1 is the x-position and c1 is the sigma of the Gaussian
                        FO.Upper        = [   100*pks(n) locs(n)+w(n)   w(n)*4]; % upper limits of the coefficients
                        FO.Lower        = [     0*pks(n) locs(n)-w(n)   w(n)/4]; % lower limits of the coefficients
                    case 'asymmG'
                        FO.StartPoint   = [      pks(n)  locs(n)      2*w(n)    0]; % start point [a1 b1 c1], where a1 is the peak amplitude, b1 is the x-position and c1 is the sigma of the Gaussian
                        FO.Upper        = [   10*pks(n)  locs(n)+w(n)   w(n)*4  1]; % upper limits of the coefficients
                        FO.Lower        = [    0*pks(n)  locs(n)-w(n)   w(n)/4 -1]; % lower limits of the coefficients
                end
                Gaussfit                = fit(xfit(:),yfit(:),fiteq,FO);
                xyPl{n}                 = [xfit(:),yfit(:)];
                %______________________________________________________________________
                % save position and amplitude of the Gaussianfit:
                Gampl(n)                = Gaussfit.a1;
                Gpos(n)                 = Gaussfit.b1;
                Gwds(n)                 = Gaussfit.c1;
                Gfits{n}                = Gaussfit; % whole fit (including uncertainty information)
                GaussfitPlot            = GaussfitPlot+feval(Gaussfit,xwin);
                %______________________________________________________________________
            end
        end
        %******************************************************************
        % Get the positions of the measured peaks above a certain threshold.
        TF0         = pks>SigThresh;
        nsig        = sum(TF0);
        SPeaksPos(jj,1:nsig+1)= [jj locs(TF0)];
        SPeaksAmp(jj,1:nsig+1)= [jj pks(TF0)];
        
        
        % Only max peak found by findpeaks. 
        [~,posmax]   = min(abs(pks-max(pks)));
        GPeaks(jj,1) = pks(posmax);
        GPeaks(jj,2) = locs(posmax);
        GPeaks(jj,3) = w(posmax);
        [~,posmax]   = min(abs(Gampl-max(Gampl)));
        GPeaks(jj,4) = Gampl(posmax);
        GPeaks(jj,5) = Gpos(posmax);
        GPeaks(jj,6) = Gwds(posmax);
        [~,posmax]   = min(abs(ywin-max(ywin)));
        ywin0 = ywin; ywin0(ywin0>thresh)=nan;
        GPeaks(jj,7) = median(maxk(ywin,3),'omitnan');
        GPeaks(jj,8) = xwin(posmax);
        % GPeaks [MaxPksAmp MaxPksPos MaxPksWdt GPksAmp GPksPos GPksWdt MaxPmeanAmp MaxPPos]
        %______________________________________________________________________
        %**************************************************************************
        %**************************************************************************
        %______________________________________________________________________
        gcf; try cla(ax1);cla(ax2);cla(ax3);cla(ax4);end % just for the plot
        %______________________________________________________________________
        % Axes 1.
        p1 = plot(ax1,xdat,ydat,'-','marker','.','color',[cm(1,:) .25],...
            'linewidth',1); hold on;
        p2 = plot(ax1,[1 1]*xwin(1),[0 max(ydat)]*1.2,'-','color',cm(2,:)); % show the window
        plot(ax1,[1 1]*xwin(end) ,[0 max(ydat)]*1.2,'-','color',cm(2,:)); % show the window
        patch(ax1,[xwin(1) xwin(1) xwin(end) xwin(end) xwin(1)],...
            [0  max(ydat)  max(ydat) 0 0]*1.2,'r','edgecolor','none','facealpha',.15)
        leg = legend(ax1,[p1 p2],'example data','considered window');
        set(leg,'interpreter','latex','box','on','location','northeast')
        
        % Axes 2.
        plot(ax2,xwin,ywin,'color',[cm(1,:) .8],'linewidth',.1); hold on;
        plot(ax2,[min(xwin) max(xwin)],[1 1]*baselineEst,'--','linewidth',1.5,'color',[cm(4,:) 1])
        title(ax2,'findpeaks Gauss fit start values','interpreter','latex')
        plot(ax2,[min(xwin) max(xwin)],[1 1]*min_height,':','color','k') % black line look at the minimum height
        if isempty(locs)
            disp('No peaks found.')
        else
            plot(ax2,([1 1].*locs(:))',(repmat([min(ywin) max(ywin)],length(locs),1))','--','color',[.5 .5 .5 .75],'linewidth',1.5 ) % peaks that were found
            plot(ax2,xwin,sum(Gauss),'color',[cm(3,:) .25],'linewidth',2) % sum of the findpeaks-Gauss, yellow curve
            for n = 1:length(locs)
                plot(ax2,xwin,Gauss(n,:),'linewidth',.1,'color',[cm(2,:) .5]);
            end
        end
        
        % Axes 3. Plot the fit (x,y) and the Gaussian fits (loop!)
        plot(ax3,xwin(:),ywin(:),'-','marker','.','linewidth',1,...
            'color',[cm(1,:) .5]); hold on % data pts considered for a fit
        plot(ax3,xwin,GaussfitPlot,'-','color',[cm(3,:) .75],'linewidth',2)
        for nn = 1:length(xyPl)
            xy = xyPl{nn};   xy1 = xyPl{1}; 
            plot(ax3,min(xy(:,1))*[1 1],[0 1.2*max(xy1(:,2))],'-','color',cm(2,:)); hold on % data pts considered for a fit
            plot(ax3,max(xy(:,1))*[1 1],[0 1.2*max(xy1(:,2))],'-','color',cm(2,:)); hold on % data pts considered for a fit
            plot(ax3,xwin,feval(Gfits{nn},xwin),'-','color',cm(2,:))
        end
        % Residuum of the Gaussfits.
        yyaxis(ax3,'right')
        try delete(py2);end
        res = ywin(:)-GaussfitPlot;%feval(Gaussfit,xfit);
        py2 = plot(ax3,xwin,res,'-','color',[cm(5,:) .5],'linewidth',.5);
        set(ax3,'ylim',[-1 1]*1.2*ceil(abs(median(maxk(res*100,20)))/100),...
            'ycolor',cm(5,:))
        ylabel(ax3,'residuum (a.u.)','interpreter','latex')
        yyaxis(ax3,'left')
        title(ax3,'Gauss fits','interpreter','latex')

        % Axes 4. Check residuum.
        histogram(ax4,res,'normalization','probability'), hold on
        histogram(ax4,ydat1,'normalization','probability')
        title(ax4,'noise and residuum of the fits','interpreter','latex')
        
        % Collect the Gaussfif values for several shots.
        plot(ax5,GPeaks(jj,2),GPeaks(jj,1),'x','markersize',8,'color',cm(1,:))
        plot(ax6,GPeaks(jj,2),GPeaks(jj,3),'x','markersize',8,'color',cm(1,:))
        
        % Axes settings.
        set([ax2 ax3],'ylim',[-10 max(ywin)*1.2],'xlim',[min(xwin) max(xwin)])
        set(ax1,'ylim',[0 max(ydat)*1.2])
        pause(.1) % otherwise can not see the control display
        %**********************************************************************
        % All peaks found in the window. Single shot info.
        Res.Gampl       = Gampl;
        Res.Gpos        = Gpos;
        Res.Gwds        = Gwds;
        Res.Gfits       = Gfits; % whole fit (including uncertainty information)
        Res.GaussfitPlot= GaussfitPlot;
        Res.xwin        = xwin;
        Res.ywin        = ywin;
        % Save all information in a data struct.
        save([wd,strrep( ['SingShot',num2str(jj),'_allPeaks_Win',num2str(round(min(xwin),2)),'-',...
            num2str(round(max(xwin),2)),'us'] ,'.',','),'.mat'],'Res')
    end
    
    
    
    
    %**********************************************************************
    if size(GPeaks,1)==1, continue, end
    GPeaks(:,9) = 1:length(GPeaks);%GPeaks = [GPeaks,(1:length(GPeaks))'];
    TF      = zeros(length(GPeaks),1);
    for ii  = 1:length(GPeaks)
        if or(sum(GPeaks(ii,1:8))==0,isnan(sum(GPeaks(ii,1:8))))
            TF(ii)=1;
        end
    end
    GPeaks(TF==1,:) = [];
    % Save main-peak information of one window for all shots into a textfile for python. 
    save([wd,strrep( ['MaxPeak_WholeMeas_Win',num2str(round(min(xwin),2)),'-',...
        num2str(round(max(xwin),2)),'us'] ,'.',','),'.txt'],'GPeaks','-ascii')
    % GPeaks [MaxPksAmp MaxPksPos MaxPksWdt GPksAmp GPksPos GPksWdt MaxPmeanAmp MaxPPos]    
    try SignalStatisticsPlot,end
    %**********************************************************************
    TFdel   = sum(SPeaksPos,1,'omitnan')==0;
    SPeaksPos(:,TFdel==1)   =[]; SPeaksAmp(:,TFdel==1)=[];
    SPeaksPos(TF==1,:)      = [];SPeaksAmp(TF==1,:)=[];
    
    SPeaks  = [SPeaksPos SPeaksAmp]; % 50 und 50
    save([wd,strrep( ['SigPeaks_WholeMeas_Win',num2str(round(min(xwin),2)),'-',...
        num2str(round(max(xwin),2)),'us'] ,'.',','),'.txt'],'SPeaks','-ascii')
    
    
end




%% Ergebnis für den Max Peak im Fenster.
TF      = zeros(length(GPeaks),1);
for ii  = 1:length(GPeaks)
    if or(sum(GPeaks(ii,:))==0,isnan(sum(GPeaks(ii,:))))
        TF(ii)=1;
    end
end
GPeaks(TF==1,:) = nan; 

figure('color','w','position',[100 100 1000 800])
ax1 = axes('position',[.1 .70 .78 .28]); hold on; box on; grid on;
ax2 = axes('position',[.1 .39 .78 .28]); hold on; box on; grid on;
ax3 = axes('position',[.1 .08 .78 .28]); hold on; box on; grid on;
set([ax1 ax2 ax3],'ticklabelinterpreter','latex','fontsize',18,'xlim',...
    [1 length(GPeaks)])
set([ax1,ax2],'xticklabel',[])

% Positions. 
p11 = plot(ax1,GPeaks(:,2),'.','markersize',12,'color',cm(1,:)); hold on
p12 = plot(ax1,GPeaks(:,5),'o','color',cm(1,:));
p13 = plot(ax1,GPeaks(:,8),'x','markersize',10,'color',cm(1,:));
leg = legend(ax1,[p11 p12 p13],'MaxPeak','Gfit','MaxMean');
set(leg,'interpreter','latex','box','on','position',[.87 .85 .1 .1]);

% Ampltidues. 
p21 = plot(ax2,GPeaks(:,1)*sqrt(2*pi).*GPeaks(:,3),'.','markersize',12,'color',cm(1,:)); hold on
p22 = plot(ax2,GPeaks(:,4)*sqrt(2*pi).*GPeaks(:,6),'o','color',cm(1,:));
p23 = plot(ax2,GPeaks(:,7),'x','markersize',10,'color',cm(1,:));
% leg = legend(ax2,[p11 p12 p13],'MaxP','Gfit','MaxMean');
% set(leg,'interpreter','latex','box','off');
plot(ax2,[0 length(GPeaks)],[1 1]*482,'-','color',cm(2,:)) % saturation

% Widths. 
p31 = plot(ax3,GPeaks(:,3),'.','markersize',12,'color',cm(1,:)); hold on
p32 = plot(ax3,GPeaks(:,6),'o','color',cm(1,:));
% leg = legend(ax3,[p11 p12],'MaxP','Gfit');
% set(leg,'interpreter','latex','box','off');

xlabel(ax3,'shot no','interpreter','latex')
ylabel(ax1,'position (m/z)','interpreter','latex')
ylabel(ax2,'amplitude (V)','interpreter','latex')
ylabel(ax3,'widths (m/z)','interpreter','latex')



h   = gcf;
fig_name = strrep([wd,'Positions_Amplitudes'],'.','');
set(h,'unit','inches');
pos = get(h,'position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
print(strcat(fig_name),'-dpng')

%% Messung 08 Proton Peak exists
TF      = cellfun(@isempty,GPeaks(:,1));
shots   = 1:size(tt,1);
p_ind   = shots(TF==0);


figure('color','w','position',[100 100 1000 600])
lagDiff = nan*ones(1,length(p_ind));
for ii = 1:length(p_ind)
    ypl = yy(p_ind(ii),:); 
    ypl(ypl>=481)=nan; % remove saturated points
    tpl = tt(p_ind(ii),:);
%     plot3(ones(length(tpl),1)*ii,tpl,ypl,'-'); hold on
    plot(tpl,ypl-ii*400,'-'); hold on
    
    % Is there a lag between the single shots?
    sref        = yy(p_ind(2),:); % hard coded spectrum no 2 as reference
    sref(sref>=480)=nan; % remove saturated points
    s1          = ypl;
    [acor,lag]  = xcorr(sref(600:1000),s1(600:1000));
    [~,I]       = max(abs(acor));
    lagDiff(ii) = lag(I);
end 

set(gca,'ticklabelinterpreter','latex','fontsize',18); box on;
% set(gca,'ylim',[.5 1])
set(gca,'xlim',[.6 1])



%% Teste asymm Gauss Funktion
% F0  = 0; % y-axis offset
    % F1  = 0; % not needed?
    % xx  = -100:100;
    % a10 = 1; % peak amplitude
    % f10 = 0.0; % Schiefheit
    % b10 = 1; % peak position
    % c10 = 10; % peak width
    % yy = F0+F1.*xx+(a10+f10.*(xx-b10)/c10)./((xx-b10).^2/c10.^2+1);
    % figure,plot(xx,yy), hold on
    % plot([xx(1) xx(end)],[0 0],'color',[.1 .1 .1 .5]) % zero-line