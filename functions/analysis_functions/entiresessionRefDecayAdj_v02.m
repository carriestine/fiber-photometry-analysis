function [FitFunction_SAVE,FitCurve_SAVE,DataDecayAdj_SAVE] = entiresessionRefDecayAdj_v02(structNames,DataRecentered,setupParam,allTrials,DecayReference)

% does not call other custom functions

%%
allTrialsID = allTrials.ID;
groupnames1 = fieldnames(allTrialsID);

groupidx = zeros(1,numel(groupnames1));

for i = 1:numel(groupnames1)
    if ~isempty(allTrialsID.(groupnames1{i}))
        groupidx(i) = 1;
    end
end

groups = groupnames1(groupidx==1);

%% Extract data for reference groups

DRefName = fieldnames(DecayReference);

for i = 1:numel(groups)
    if ismember(groups{i},DRefName)
        DecayRef.(groups{i}) = cell2mat(DecayReference.(groups{i}));
    end
end


exptGroups = fieldnames(DecayRef);
RefNames = cell(1,numel(exptGroups));
for i = 1:numel(exptGroups)
    RefNames{i} = DecayRef.(groups{i});
end

for t = 2:numel(structNames)
    IsolatedRef.(structNames{t}) = cell(1,numel(RefNames));
    for i = 1:numel(RefNames)
        refidx = ismember(allTrials.Grouped,allTrialsID.(RefNames{i}));
        IsolatedRef.(structNames{t}){i} = DataRecentered.(structNames{t})(:,refidx);
    end
end


%% Average data in reference groups

for t = 2:numel(structNames)
    avgRef.(structNames{t}) = zeros(size(DataRecentered.(structNames{t}),1),numel(IsolatedRef.(structNames{t})));
    for i = 1:numel(IsolatedRef.(structNames{t}))
        avgRef.(structNames{t})(:,i) = mean(IsolatedRef.(structNames{t}){i},2);
    end
end
%% Fit curve to average of each reference group, store as columns in an array

tempX = DataRecentered.Dts(:,1);
for t = 2:numel(structNames)
   % FitCurve1.(structNames{t}).(RefNames{i}) = zeros(size(DataRecentered.(structNames{t}),1),numel(IsolatedRef.(structNames{t})));
    for i = 1:numel(RefNames)
        FitCurve1.(structNames{t}).(RefNames{i}) = zeros(size(DataRecentered.(structNames{t}),1),1);
        if contains(setupParam.CurveType,'poly')
            FitFunction1.(structNames{t}).(RefNames{i}) = polyfit(tempX,avgRef.(structNames{t})(:,i),4);
            FitCurve1.(structNames{t}).(RefNames{i}) = polyval(FitFunction1.(structNames{t}).(RefNames{i}),tempX);
        elseif contains(setupParam.CurveType,'exp')
            FitFunction1.(structNames{t}).(RefNames{i}) = fit(tempX,avgRef.(structNames{t})(:,i),'exp2');
            FitCurve1.(structNames{t}).(RefNames{i}) = FitFunction1.(structNames{t}).(RefNames{i})(tempX);
        end
    end
end

%% Create FitCurve matrix where each column corresponds to the fit curve that should be subtracted from corresponding column of data matrix
for t = 2:numel(structNames)
    FitCurveCompiled.(structNames{t}) = zeros(size(DataRecentered.Dts));
    for i = 1:numel(exptGroups)
        colidx = ismember(allTrials.Grouped,allTrialsID.(exptGroups{i}));
        for p = find(colidx==1)
            FitCurveCompiled.(structNames{t})(:,p) = FitCurve1.(structNames{t}).(DecayRef.(exptGroups{i}));
        end
    end
    
    for r = 1:numel(RefNames)
        colidx = ismember(allTrials.Grouped,allTrialsID.(RefNames{r}));
        for p = find(colidx==1)
        	FitCurveCompiled.(structNames{t})(:,p) = FitCurve1.(structNames{t}).(RefNames{r});
        end
    end
end
    
%%  Subtract fit curve from recentered data

DataDecayAdj.signal = DataRecentered.signal - FitCurveCompiled.signal;
DataDecayAdj.reference = DataRecentered.reference - FitCurveCompiled.reference;
DataDecayAdj.sigsubref = DataDecayAdj.signal - DataDecayAdj.reference;


%%
FitFunction_SAVE = FitFunction1;
FitCurve_SAVE = FitCurveCompiled;
DataDecayAdj_SAVE = DataDecayAdj;






