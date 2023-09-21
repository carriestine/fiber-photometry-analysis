function eventwindowPhotomPlot_SigRef_v02(x,y,eb,setupParam,TrialParam,allTrials,ProcessedEventWindowData)

% calls 1 other custom function: photomplotTitle_v_02
%%
allTrialsID = allTrials.ID;

if setupParam.CompareGroups == 1
    GroupParam = TrialParam.GroupParam;
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
%     refGroups = refGroups1(ismember(refGroups1,groups));
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

% VN1 = varnamesDATA(~contains(varnamesDATA,'Z'));
% VN = VN1(~contains(VN1,'dFF'));

VN1 = {'signal' 'reference' 'sigsubref'};
VN2 = {'Zbytrial_sig' 'Zbytrial_ref' 'Zbytrial_sigsubref'};
VNall = {VN1 VN2}; 

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

% Setup variable names/color
legID = {'signal' 'reference' 'sig - ref'};
choosecolor = {'#0c9cf8' '#330672' 'g'};

tempX = x.Dts;

if setupParam.CompareGroups == 0 || setupParam.CombineGroups == 1
    for a = 1:numel(VNall)
        VN = VNall{a};
        figure(10+a)
        numcol = numel(ST)*3;
        for t = 1:numel(ST)
            spSpot = [1,2,3]+(3*(t-1));
            subplot(2,numcol,spSpot)

            for i = 1:numel(VN)
                [currenttitle] = photomplotTitle_v02(ST(t),VN(i));

                tempY = y.all.(ST{t}).(VN{i});
                tempEBlo = eb.all.lo.(ST{t}).(VN{i});
                tempEBhi = eb.all.hi.(ST{t}).(VN{i});

                leg(i) = plot(tempX,tempY,'Color',string(choosecolor(i)),'LineWidth',2); %#ok<AGROW>
                hold on
                    filledEBtemp = patch([tempX fliplr(tempX)],[tempEBlo fliplr(tempEBhi)],'k');
                    alpha(filledEBtemp,0.2)
                    filledEBtemp.EdgeAlpha = 0;
                    set(filledEBtemp,'edgecolor',string(choosecolor(i)),'facecolor',string(choosecolor(i)))
            end
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
            title(sprintf('Combined %s %s sig and ref',currenttitle.ST, currenttitle.Calc))
            legend(leg,legID)
            hold off

            % plot heat map
            for i = 1:numel(VN)
                hmSpot = (numcol+i) + (3*(t-1));
                subplot(2,numcol,hmSpot)

                [currenttitle] = photomplotTitle_v02(ST(t),VN(i));

                tempHM = transpose(ProcessedEventWindowData.(ST{t}).(VN{i})(:,expidx1:expidx2));
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
                
        end

        sname = 'CarrieCustom';
        S = hgexport('readstyle',sname);
    %     style.Format = 'png';
    %     style.ApplyStyle = '1';
        hgexport(gcf,'file_name',S,'applystyle',true);

        hold off
    end
end


%% Plot by group
if setupParam.CompareGroups == 1
    for g = 1:numel(plotgroups)
        for a = 1:numel(VNall)
            VN = VNall{a};
            figure(10*(g+1) + a)

            group_expidx1 = GroupParam.(plotgroups{g}).eventindex(1,1);
            group_expidx2 = GroupParam.(plotgroups{g}).eventindex(2,end);

            numcol = numel(ST)*3;
            for t = 1:numel(ST)
                spSpot = [1,2,3]+(3*(t-1));
                subplot(2,numcol,spSpot)

                for i = 1:numel(VN)
                    [currenttitle] = photomplotTitle_v02(ST(t),VN(i));

                    tempY = y.(plotgroups{g}).(ST{t}).(VN{i});
                    tempEBlo = eb.(plotgroups{g}).lo.(ST{t}).(VN{i});
                    tempEBhi = eb.(plotgroups{g}).hi.(ST{t}).(VN{i});

                    legGroup(i) = plot(tempX,tempY,'Color',string(choosecolor(i)),'LineWidth',2); %#ok<AGROW>
                    hold on
                        filledEBtemp = patch([tempX fliplr(tempX)],[tempEBlo fliplr(tempEBhi)],'k');
                        alpha(filledEBtemp,0.2)
                        filledEBtemp.EdgeAlpha = 0;
                        set(filledEBtemp,'edgecolor',string(choosecolor(i)),'facecolor',string(choosecolor(i)))
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
                legend(legGroup,legID)
                hold off

                % plot heat map
                for i = 1:numel(VN)
                    hmSpot = (numcol+i) + (3*(t-1));
                    subplot(2,numcol,hmSpot)

                    [currenttitle] = photomplotTitle_v02(ST(t),VN(i));

                    tempHM = transpose(ProcessedEventWindowData.(ST{t}).(VN{i})(:,group_expidx1:group_expidx2));
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
                            ticklab{s} = append('Animal ',num2str(ticklab1(s))); 
                        end

                        xlim(xl)
                        set(gca,'ytick',tickpts,'yticklabel',ticklab)
                        title(sprintf('Heat Map: %s %s',currenttitle.ST,currenttitle.Data))
                end

            end

            sname = 'CarrieCustom';
            S = hgexport('readstyle',sname);
        %     style.Format = 'png';
        %     style.ApplyStyle = '1';
            hgexport(gcf,'file_name',S,'applystyle',true);

            hold off
        end
    end
end


%% Scale all axes across groups to be the same for simple comparison
 for a = 1:numel(VNall)
    for g = 1:numel(plotgroups)
        figure(10*(g+1) + a)
        for t = 1:numel(ST)
            spSpot = [1,2,3]+(3*(t-1));
            subplot(2,numcol,spSpot)
            YLgroup.(ST{t})(g,1:2) = ylim;
        end
    end
 end

for g = 1:numel(plotgroups)
    for t = 1:numel(ST)
        YLuse.(ST{t})(1,1) = min(YLgroup.(ST{t})(:,1));
        YLuse.(ST{t})(1,2) = max(YLgroup.(ST{t})(:,2));
    end
end

for a = 1:numel(VNall)
    for g = 1:numel(plotgroups)
        figure(10*(g+1) + a)
        for t = 1:numel(ST)    
            spSpot = [1,2,3]+(3*(t-1));
            subplot(2,numcol,spSpot)
            ylim(YLuse.(ST{t}))
        end
    end
end











