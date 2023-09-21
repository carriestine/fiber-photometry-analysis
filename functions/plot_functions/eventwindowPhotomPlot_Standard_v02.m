function eventwindowPhotomPlot_Standard_v02(x,y,eb,setupParam,TrialParam,allTrials,ProcessedEventWindowData,varnamesDATA)

% calls 2 other custom functions: photomplotTitle_v02, photomplotColor_v02

%%
allTrialsID = allTrials.ID;

if setupParam.CompareGroups == 1
    groupnames1 = fieldnames(allTrialsID);
    groupidx = zeros(1,numel(groupnames1));
        for i = 1:numel(groupnames1)
            if ~isempty(allTrialsID.(groupnames1{i}))
                groupidx(i) = 1;
            end
        end
    groups = groupnames1(groupidx==1);
    
    expGroups1 =  ['A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'];
    expGroups = expGroups1(ismember(expGroups1,groups));
    
    numexpAnimals1 = zeros(1,numel(expGroups));
    for d = 1:numel(expGroups)
        numexpAnimals1(d) = numel(allTrialsID.(expGroups(d)));
    end
    numexpAnimals = sum(numexpAnimals1);
        expidx1 = TrialParam.eventindex(1,1);
        expidx2 = TrialParam.eventindex(2,numexpAnimals);
        anidx1 = TrialParam.animalindex(1);
        anidx2 = TrialParam.animalindex(numexpAnimals);
    
    refGroups1 = ['T'; 'U'; 'V'; 'W'; 'X'; 'Y'; 'Z'];
    plotgroups = groups(~ismember(groups,refGroups1));
    
    legendID1 = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
    legID = legendID1(find(ismember(expGroups1,plotgroups))); %#ok<FNDSB>
    numgroups = numel(plotgroups);
%     refGroups = refGroups1(ismember(refGroups1,groups));
else
    numgroups = 0;
end

%% Setup data name/timing parameters
All = y.all;
%ST1 = {'Raw'; 'Recentered'; 'DecayAdj'};
ST1 = {'DecayAdj'};
inUseST = fieldnames(All);
STidx = zeros(1,numel(ST1));
    for i = 1:numel(ST1)
        if ismember(ST1(i),inUseST)
            STidx(i) = 1;
        end
    end
ST = ST1(STidx==1);

VNnone = {'signal' 'sigsubref'};

VN1 = {'dFFentire_sig' 'dFFentire_sigsubref' 'Zentire_sig' 'Zentire_sigsubref'};
VNentire = varnamesDATA(ismember(varnamesDATA,VN1));


VN2 = {'dFFbytrial_sig' 'dFFbytrial_sigsubref' 'Zbytrial_sig' 'Zbytrial_sigsubref'};
VNbytrial = varnamesDATA(ismember(varnamesDATA,VN2));

VNall = {VNnone VNentire VNbytrial};

eventtime = 0;
eventduration = TrialParam.eventduration(1);
timewindow = setupParam.timewindow;

% setup x limits
xl = [0 0];
if setupParam.viewtimewindow == 1
    xl(1) = -setupParam.viewtimewindowVal;
    xl(2) = setupParam.viewtimewindowVal;
else
    xl(1) = -timewindow;
    xl(2) = timewindow;
end

%% Plot combined signal, reference, and signal - reference on one plot

[choosecolor] = photomplotColor_v02(setupParam,numgroups);

tempX = x.Dts;

if setupParam.CompareGroups == 0 || setupParam.CombineGroups == 1
   for v = 1:numel(VNall)
       figure(20 + v)
       VN = VNall{v};
 
       numcol = numel(VN);
       numrow = 2;
       for i = 1:numel(VN)
          spSpot = i;
          subplot(numrow,numcol,spSpot)

          [currenttitle] = photomplotTitle_v02(ST(1),VN(i));

          tempY = y.all.DecayAdj.(VN{i});
          tempEBlo = eb.all.lo.DecayAdj.(VN{i});
          tempEBhi = eb.all.hi.DecayAdj.(VN{i});

          plot(tempX,tempY,'Color',string(choosecolor(2)),'LineWidth',2)
          hold on
               filledEBtemp = patch([tempX fliplr(tempX)],[tempEBlo fliplr(tempEBhi)],'k');
                alpha(filledEBtemp,0.2)
                filledEBtemp.EdgeAlpha = 0;
                set(filledEBtemp,'edgecolor',string(choosecolor(2)),'facecolor',string(choosecolor(2)))
                
        xlim(xl)
        yl = ylim;
        tempYL = ((floor(yl(1)) - 1):(round(yl(2)) + 1));
            if any(ismember((tempYL),0))
                line([tempX(1) tempX(end)], [0 0],'Color','k');
            end

        line([eventtime eventtime],[yl(1) yl(2)],'Color','k');
        line([(eventtime + eventduration) (eventtime + eventduration)],[yl(1) yl(2)],'Color','k');
        ylim(yl)

        ylabel(sprintf('%s',currenttitle.Yaxis));
        xlabel('Time (s)')
        title(sprintf('Combined %s %s %s',currenttitle.ST, currenttitle.Calc,currenttitle.Data))
       end

        

        % plot heat map
                for i = 1:numel(VN)
                    hmSpot = (numcol + i);
                    subplot(numrow,numcol,hmSpot)

                    [currenttitle] = photomplotTitle_v02(ST(1),VN(i));

                    tempHM = transpose(ProcessedEventWindowData.DecayAdj.(VN{i})(:,expidx1:expidx2));
                        imagesc(tempX,1:size(tempHM,1),tempHM,yl)
                        hold on
                        cb = colorbar;
                        title(cb,currenttitle.Yaxis)
                        L2 = line([0 0], [0 size(tempHM,1)+1]);
                        set(L2,'Color','white','LineWidth',3)

                        for p = 1:numexpAnimals
                            hmYlim = TrialParam.eventindex(2,p);
                            L = line([-timewindow timewindow], [hmYlim+0.5 hmYlim+0.5]);
                            set(L,'Color','white','LineWidth',3)
                        end

                        tickpts1 = TrialParam.eventindex(:,anidx1:anidx2);
                        tickpts = zeros(1,size(tickpts1,2));
                        ticklab1 = zeros(1,size(tickpts1,2));
                        for s = 1:size(tickpts1,2)
                            tickpts(s) = tickpts1(2,s) - ((tickpts1(2,s) - tickpts1(1,s))./2);
                            ticklab1(s) = TrialParam.animalindex(s);
                            ticklab{s} = append('Animal ',num2str(ticklab1(s))); %#ok<AGROW>
                        end

                    xlim(xl)
                    set(gca,'ytick',tickpts,'yticklabel',ticklab)
                    title(sprintf('Heat Map: %s %s',currenttitle.ST,currenttitle.Data))
                end

            sname = 'CarrieCustom';
            S = hgexport('readstyle',sname);
            hgexport(gcf,'file_name',S,'applystyle',true);

            hold off
   end
end

%% Plot by group

if setupParam.CompareGroups == 1
    for v = 1:numel(VNall)
        figure(30 + v)
        VN = VNall{v};
        
        numcol = numel(VN);
        numrow = 1;
        
        for i = 1:numel(VN)
            spSpot = i;
            subplot(numrow,numcol,spSpot)
            
            [currenttitle] = photomplotTitle_v02(ST(1),VN(i));
            leg = zeros(1,numel(plotgroups));
            for g = 1:numel(plotgroups)
                tempY = y.(plotgroups{g}).DecayAdj.(VN{i});
                tempEBlo = eb.(plotgroups{g}).lo.DecayAdj.(VN{i});
                tempEBhi = eb.(plotgroups{g}).hi.DecayAdj.(VN{i});
                
                leg(g) = plot(tempX,tempY,'Color',string(choosecolor(g)),'LineWidth',2); 
                hold on
                filledEBtemp = patch([tempX fliplr(tempX)],[tempEBlo fliplr(tempEBhi)],'k');
                alpha(filledEBtemp,0.2)
                filledEBtemp.EdgeAlpha = 0;
                set(filledEBtemp,'edgecolor',string(choosecolor(g)),'facecolor',string(choosecolor(g)))
            end
            
            xlim(xl)
            yl = ylim;
            tempYL = ((floor(yl(1)) - 1):(round(yl(2)) + 1));
            if any(ismember((tempYL),0))
                line([tempX(1) tempX(end)], [0 0],'Color','k');
            end
            
            line([eventtime eventtime],[-1000 1000],'Color','k');
            line([(eventtime + eventduration) (eventtime + eventduration)],[-1000 1000],'Color','k');
            ylim(yl)
            
            ylabel(sprintf('%s',currenttitle.Yaxis));
            xlabel('Time (s)')
            title(sprintf('Group %s: %s %s sig and ref',plotgroups{g},currenttitle.ST,currenttitle.Calc))
            legend(leg,legID)
        end
        
        sname = 'CarrieCustom';
        S = hgexport('readstyle',sname);
        %     style.Format = 'png';
        %     style.ApplyStyle = '1';
        hgexport(gcf,'file_name',S,'applystyle',true);
        
        hold off
    end
end

 
