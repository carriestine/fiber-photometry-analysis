%% IN PROGRESS CODE
% This version fits a curve to the recentered traces to account for
% bleaching decay prior to performing Z-score and DF/F calculations. (Decay
% adjustment prior to normalization). It is also updated to clean up the 
% miscellaneous 'scale reference to signal' processing that was not
% fully functional, as well as reorganizing the order of variables in the 
% adjust section.

%% Reset MatLab workspace - clears all variables and the command window

clear variables;  % clear all variables
close all;  % close all open graphs
clc   % clear command window

%% EDIT THESE FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enter trial info
spreadsheetName = 'nameGenerator_NOPLightGCaMP_v02.csv';

T = [] - 1; %ADJUST input rows of nameGenerator spreadsheet trials that correspond to your first reference group
U = [] - 1; %ADJUST input rows of nameGenerator spreadsheet trials that correspond to your second reference group
V = [] - 1; %ADJUST see above
W = [] - 1; %ADJUST see above
X = [] - 1; %ADJUST see above
Y = [] - 1; %ADJUST see above
Z = [] - 1; %ADJUST see above

A = [2:4] - 1; %#ok<*NBRAK> %ADJUST row numbers for trials averaged in group 1
DecayReference.A = {'T'}; %ADJUST change to whichever reference group corresponds to experimental group 1

B = [23:25] - 1; %ADJUST row numbers for trials averaged in group 2
DecayReference.B = {'T'}; %ADJUST change to whichever reference group corresponds to experimental group 2

C = [] - 1; %ADJUST see above
DecayReference.C = {'T'}; %ADJUST see above

D = [] - 1; %ADJUST see above
DecayReference.D = {'W'}; %ADJUST see above

E = [] - 1; %ADJUST see above
DecayReference.E = {'T'}; %ADJUST see above

F = [] - 1; %ADJUST see above
DecayReference.F = {'Y'}; %ADJUST see above

G = [] - 1; %ADJUST see above
DecayReference.G = {'Z'}; %ADJUST see above

H = [] - 1; %ADJUST see above
DecayReference.H = {'Z'}; %ADJUST see above

% Change signal and isosbestic channel identifier
    
    sigID = {'raw470'};
    
    refID = {'raw405'};

% Set group/averaging/smoothing parameters

    RefDecayAdjust = 0; % ADJUST to 1 if you want to match exp trials to reference trials based on curve clustering as the method for decay adjustment
    
    CompareGroups = 1; % ADJUST to 1 if you are comparing multiple groups (A vs B, etc), 0 if you are just looking at trials in A
    CombineGroups = 0; % ADJUST to 1 if you want to look at the combination of all groups post decay adjustment
    
    averageall = 1; % ADJUST to 1 if you want to average all trials for all animals, or 0 if you want to average trials within an animal then average those averages
    
    resamplefactor = 300; % ADJUST to factor you want to resample data at
    
    CurveType = 'poly'; % ADJUST to 'exp' or 'poly' depending on if you want to fit exponential or polynomial curves for decay adjustment
    
    ScaleRef = 1; % ADJUST to 1 to detrend bleaching decay by fitting (scaling) the reference to the signal
       ScaleRef_PostDecay = 0; % ADJUST to 1 to fit (scale) the reference to the signal AFTER fitting/subtracting a curve to account for decay
        scaleref_BLper = 1; % ADJUST to 1 to scale ref to sig using JUST the BL period bleaching. Good for entire session decay with pharmacology, but for transient events set this to 0 to scale using entire trace

    doSmooth = 1; % ADJUST to 1 if want to apply smoothing to the trace
    smooth_win = 1; % ADJUST period of time for smoothing, in seconds


% Set timing parameters

    manualsessionlength = 0; % ADJUST if sessions compared are not the same length
    newsessionlength = 1800;% ADJUST to lowest common session length, in seconds
    
    BLperiod = [21 570]; %ADJUST to period of time where you want to baseline the trace to (in seconds, ie 1 to 260s when injection is at 300s)
    
    starttime = 20; % ADJUST to time (in seconds) that you want to treat as the beginning of the session
    
    %Event window specific stuff below, doesn't need to be adjusted for
    %this entire session script (just keeping it here to make updating
    %either version easier)
    timewindow = 60; %ADJUST +/- window (seconds) around each event that you want to pull out data for (must contain idxBase values in range)
    
    viewtimewindow = 0; % ADJUST if you want to plot a smaller subset of time than the entire time window of data you'll be analyzing
    viewtimewindowVal = 20; % ADJUST +/- window (seconds) that you want to see in the final plot
    
    manualtimestamps = 0; % ADJUST if you want to override the session's timestamps to enter uniform timestamps manually for all sessions
    newtimestamps = [120;240;360;480]; %ADJUST to array of timestamps that you want to use instead

    idxBase = [-30 -15]; % ADJUST to window of time (seconds) relative to each event that you want to use as each event's baseline window 

% Set plot color parameters

    comparisonOptions = {'1: green gradient' '2: blue gradient' '3: red gradient' '4: gray gradient' '5: paired multicolor' '6: 405/470 vs 435/490' '7: choose colors in photomplotColor'};
    comparisonType = 1; %ADJUST to number corresponding to comparisonOption you want the colors to be. Gradients are light to dark

    PlotIndividual = 0; %ADJUST if you want to see graphs for each ind animal
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Establish paths, setup spreadsheet info as cell array

path_to_functions = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/fiber-photometry-analysis/functions';
 addpath(path_to_functions);
 
path_to_functions2 = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/fiber-photometry-analysis/functions/analysis_functions';
 addpath(path_to_functions2);

path_to_functions3 = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/fiber-photometry-analysis/functions/plot_functions';
 addpath(path_to_functions3); 
 
path_to_spreadsheet = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/fiberphotometry_DATA/name_spreadsheets';
 addpath(path_to_spreadsheet);

    nameGenerator = readtable(spreadsheetName); 
    nameGen = table2cell(nameGenerator); % turn spreadsheet with trial info into a cell array
 
%% Prep variables in format needed for functions
    clear allTrials; clear setupParam;
    allTrials.Grouped = [A B C D E F G H T U V W X Y Z];
    allTrials.ID.A = A; allTrials.ID.B = B; allTrials.ID.C = C; allTrials.ID.D = D; allTrials.ID.E = E; allTrials.ID.F = F; allTrials.ID.G = G; allTrials.ID.H = H; allTrials.ID.T = T; allTrials.ID.U = U; allTrials.ID.V = V; allTrials.ID.W = W; allTrials.ID.X = X; allTrials.ID.Y = Y; allTrials.ID.Z = Z; allTrials.ID.R = [T U V W X Y Z];
    
   
    setupParam_NAMES = {'sigID' 'refID' 'CompareGroups', 'CombineGroups', 'averageall', 'resamplefactor', 'CurveType', 'RefDecayAdjust', 'ScaleRef', 'ScaleRef_PostDecay', 'scaleref_BLper', 'doSmooth', 'smooth_win', 'manualsessionlength', 'newsessionlength', 'manualtimestamps', 'newtimestamps', 'starttime', 'timewindow', 'viewtimewindow', 'viewtimewindowVal', 'BLperiod', 'idxBase', 'comparisonType', 'PlotIndividual'};
    setupParam_VALUES = {sigID refID CompareGroups CombineGroups averageall resamplefactor CurveType RefDecayAdjust ScaleRef ScaleRef_PostDecay scaleref_BLper doSmooth smooth_win manualsessionlength newsessionlength manualtimestamps newtimestamps starttime timewindow viewtimewindow viewtimewindowVal BLperiod idxBase comparisonType PlotIndividual};

    setupParam = cell2struct(setupParam_VALUES,setupParam_NAMES,2);

    
    
    
    
    
    
    
%% Combine time series info for all trials (together and by group) into a single struct

[PreprocessedData,TrialParam,DataVersions,varnames,Fs] = entiresessionStackTrials_v02(nameGen,allTrials,setupParam,DecayReference); % original was entiresessionStackTrialsV5

varnamesDATA = varnames(2:end);

%% Plot individual traces
if setupParam.PlotIndividual == 1
    plotIndividualEntire_v02(PreprocessedData,allTrials,setupParam,TrialParam.FitCurve) % original was plotIndividualEntire (no version)
end
    
%% Scale isosbestic to signal using LLS after decay adjustment
if ScaleRef_PostDecay == 1
    [scaledRef,scaledSigSubRef] = scaleRefToSig_postDecay_v02(PreprocessedData.DecayAdj.Dts,PreprocessedData.DecayAdj,setupParam,TrialParam.delay,allTrials); % original was scaleRefToSig_postDecayV4
        PreprocessedData.DecayAdj.reference = scaledRef;
        PreprocessedData.DecayAdj.sigsubref = scaledSigSubRef;
end

%% Calculate ∆F/F of entire session using fixed BL period in the session

clc
[ProcessedData,varnamesDATA,Dts] = fixedBLperDFF_v02(PreprocessedData,TrialParam,varnames,varnamesDATA); % original was fixedBLperDFFV4


%% Calculate Z-score of entire session using fixed BL period in the session

[ProcessedData,varnamesDATA] = fixedBLperZscore_v02(ProcessedData,Dts,TrialParam,varnames,varnamesDATA); % original was fixedBLperZscoreV4



%% Generate within group averages and calculate SEM

[DataAVG] = entiresessionAverage_v02(ProcessedData,varnamesDATA,setupParam,TrialParam,allTrials); % original was entiresessionAverageV2

%% Format final processed data for plotting

[x,y,eb] = entiresessionPlotFormat_v02(DataAVG,Dts,varnamesDATA,setupParam,TrialParam,allTrials); % original was entiresessionPlotFormatV2

%% Obtain average signal for each animal in 5 min time bins

BinParam.type  = {'DecayAdj'}; % can be: 'Raw' 'Recentered' 'DecayAdj'
BinParam.process = 'Zentire'; % can be: 'none' 'dFFentire' 'Zentire'
BinParam.data = 'sigsubref'; % can be: 'sig', 'sigsubref'

[BIN_Table,binCheck] = entiresessionAvg5minBins_v02(ProcessedData,setupParam,BinParam,TrialParam,x); % original was Avg5minBins

%% Obtain area under curve for each animal

AUCParam.type = {'DecayAdj'}; % can be: 'Raw' 'Recentered' 'DecayAdj'
AUCParam.process = 'dFFentire'; % can be: 'none' 'dFFentire' 'Zentire'
AUCParam.data = 'sigsubref'; % can be: 'sig', 'sigsubref'

[AUC_Table] = entiresessionAUC_v02(ProcessedData,setupParam,AUCParam,TrialParam,x); % original was entiresessionAUC
%% Plot signal and reference together

%entiresessionPhotomPlot_SigRef_v02(x,y,eb,setupParam,TrialParam,allTrials); % original was entiresessionPhotomPlot_SigRefV2



%% Plot standard figures together

%entiresessionPhotomPlot_Standard_v02(x,y,eb,setupParam,TrialParam,allTrials); % original was entiresessionPhotomPlot_StandardV2

%% Plot figures as you choose
loneFigure.fignum = 101;
loneFigure.type = {'DecayAdj'}; % can be: 'Raw' 'Recentered' 'DecayAdj'
loneFigure.data = 'Zentire'; % can be: 'none' 'dFFentire' 'Zentire'
entiresessionPhotomPlot_LoneFig_v02(x,y,eb,setupParam,TrialParam,allTrials,loneFigure); % original was entiresessionPhotomPlot_LoneFigV2



%% Plot figures as you choose
loneFigure.fignum = 201;
loneFigure.type = {'Raw'}; % can be: 'Raw' 'Recentered' 'DecayAdj'
loneFigure.data = 'none'; % can be: 'none' 'dFFentire' 'Zentire'
entiresessionPhotomPlot_LoneFig_v02(x,y,eb,setupParam,TrialParam,allTrials,loneFigure);














%% SAVE CUSTOM FIGURE FOR PRESENTATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savePLOT = 'Z'; % can be: 'Z' 'Raw'
figname = 'NOPLightGqDREADD_470-405_Ro10mgkg_ScaledRef_none_159432.eps';
figtitle = 'NOPLight Gq DREADD Ro 10 mg/kg (470-405)';
linecolor = '#b093e1'; % #b093e1 (purple)     #6180a5 (gray)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tempX = x.Dts;

if contains(savePLOT,'Z')
    tempY = y.A.Recentered.reference;
    tempEBlo = eb.A.lo.Recentered.reference;
    tempEBhi = eb.A.hi.Recentered.reference;
    yl = [-30 20];
elseif contains (savePLOT,'dFF')
    tempY = y.all.DecayAdj.dFFentire_sigsubref;
    tempEBlo = eb.all.lo.DecayAdj.dFFentire_sigsubref;
    tempEBhi = eb.all.hi.DecayAdj.dFFentire_sigsubref;
    yl = [-50 50];
elseif contains(savePLOT,'Raw')
    tempY = y.all.DecayAdj.sigsubref;
    tempEBlo = eb.all.lo.DecayAdj.sigsubref;
    tempEBhi = eb.all.hi.DecayAdj.sigsubref;
    yl = [-1 1];
end

eventtime = TrialParam.rawtimestamps{1}(1);
eventduration = TrialParam.eventduration(1);


%% Plot raw signal/isosbestic on left and (signal - isosbestic) on right, colored by group
groupid = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
numgroups = size(TrialParam.GroupParam.allgroups,1);
[choosecolor] = photomplotColor_v02(setupParam,numgroups);

for k = 1:numgroups
    for g = TrialParam.GroupParam.(groupid{k}).animalindex(1):TrialParam.GroupParam.(groupid{k}).animalindex(end)
        figure(1000+g)
        subplot(1,2,1)
            plot(Dts(:,g),ProcessedData.Recentered.signal(:,g),'c')
            hold on
            plot(Dts(:,g),ProcessedData.Recentered.reference(:,g), 'Color', '#b093e1')
                yl = ylim;
                line([600 600], [yl(1) yl(2)],'Color','k');
                line([2400 2400], [yl(1) yl(2)],'Color','k');
                ylim(yl);
                
        subplot(1,2,2)
            plot(Dts(:,g),ProcessedData.Recentered.sigsubref(:,g), 'Color', string(choosecolor(k)))
                yl = ylim;
                line([600 600], [yl(1) yl(2)],'Color','k');
                line([2400 2400], [yl(1) yl(2)],'Color','k');
                ylim(yl);
                
    sname = 'CarrieCustom';
    S = hgexport('readstyle',sname);

    hgexport(gcf,'file_name',S,'applystyle',true);
    end
end



