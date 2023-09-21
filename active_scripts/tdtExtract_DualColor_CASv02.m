% Code written by Carrie Stine, parts adapted from code written by Christian Pedersen
% Michael Bruchas Lab - UW

% Extracts photometry data from TDT Data Tank and outputs time series as
% .mat file

% REQUIRES tdt2mat.m file to be in same folder (filepath)

%% Reset MatLab workspace - clears all variables and the command window

clear all;  % clear all variables
close all;  % close all open graphs
clc   % clear command window



%%
A = [125];
for a1 = A %change this range to rows of your nameGenerator file you want to run MINUS ONE (table doesn't import headers)
   %%%%%%%%%%%%%%%%%%%%%%  EDIT THESE FIELDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    % point to tanks (enter file path)
    %path_to_data = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository'; %change to point to where your data is
    path_to_data = '/Users/castine/Dropbox (Bruchas Lab)/Carrie Data/Bruchas Lab/pnVTA nociceptin/DATA/NOPLight validation/Tanks';
    % set up Tank name, variables to extract
    tankdir = [path_to_data];
    
    path_to_data2 = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/fiber-photometry-analysis';
    addpath(path_to_data2);

    path_to_data3 = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/fiber-photometry-analysis/functions';
    addpath(path_to_data3);

    path_to_data4 = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/fiber-photometry-analysis/functions/extract_functions';
    addpath(path_to_data4);

    path_to_spreadsheet = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/fiberphotometry_DATA/name_spreadsheets';
    addpath(path_to_spreadsheet);
    
    
    nameGenerator = readtable('nameGenerator_NOPLightGCaMP.csv'); %enter name of your excel sheet containing naming info
    nameGen = table2cell(nameGenerator);

    plotFigures = 0; %change to 1 if you want to see the plots
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    a = a1-1;
    
    tankname = nameGen{a,1};
    
    blockname = nameGen{a,2};
    
    tempID = string(nameGen(a,:));
    
    filename1 = strjoin(tempID(3:5),'');
    filename2 = append(filename1,'.',"mat");
    filename = char(filename2);

    figname1 = append(filename1,'.',"eps");
    figname = char(figname1);

    TTLname1 = append(filename1,"_timestamps",'.',"txt");
    TTLname = char(TTLname1);

%%
%if(a,10)

FPstatic = 'Static';
FPstaticTEMP = 'StaticTEMP'; %temporary storeIDs while photodetector A is broken
FPNAPE = 'NAPE';
FPTwoMice1 = 'StaticTWOMICE_FP1';
FPTwoMice2 = 'StaticTWOMICE_FP2';
FPTwoMice1TEMP = 'StaticFP1TEMP';
FPTwoMice2TEMP = 'StaticFP2TEMP';
BruchasCartDual470 = 'BruchasCartDual470';
BruchasCartDual490 = 'BruchasCartDual490';
TTL = 'TTL';

if strcmp(nameGen{a,10},FPstatic) 
    storenames = {'465A'}; % name of stores to extract from TDT (usu. 4-letter code)
    %LMag is the demodulated data, may also have other timestamps etc
    storenames2 = {'405A'};
    storenames3 = {'560B'};
    
elseif strcmp(nameGen{a,10},FPstaticTEMP)
    storenames = {'465C'};
    storenames2 = {'405C'};
    storenames3 = {'560B'};

elseif strcmp(nameGen{a,10},FPNAPE)
    storenames = {'470A'};
    %storenames2 = {'405A'};
    storenames2 = {'560B'};
    storenames3 = {'560B'};
elseif strcmp(nameGen{a,10},FPTwoMice1)
    storenames = {'405A'};
    storenames2 = {'465A'};
    storenames3 = {'465A'};
elseif strcmp(nameGen{a,10},FPTwoMice2)
    storenames = {'405C'};
    storenames2 = {'465C'};
    storenames3 = {'465C'};
elseif strcmp(nameGen{a,10},FPTwoMice1TEMP)
    storenames = {'405B'};
    storenames2 = {'465B'};
    storenames3 = {'465B'};
elseif strcmp(nameGen{a,10},FPTwoMice2TEMP)
    storenames = {'405C'};
    storenames2 = {'465C'};
    storenames3 = {'465C'};
elseif strcmp(nameGen{a,10},BruchasCartDual470)
    storenames = {'470A'};
    storenames2 = {'405A'};
    storenames3 = {'565B'};
elseif contains(nameGen{a,10},BruchasCartDual490)
    storenames = {'490A'};
    storenames2 = {'435A'};
    storenames3 = {'565B'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract
for k = 1:numel(storenames)
  storename = storenames{k};
  S{k} = tdt2mat(tankdir, tankname, blockname, storename);
end

for k = 1:numel(storenames2)
  storename2 = storenames2{k};
  S2{k} = tdt2mat(tankdir, tankname, blockname, storename2);
end


for k = 1:numel(storenames3)
  storename3 = storenames3{k};
  S3{k} = tdt2mat(tankdir, tankname, blockname, storename3);
end

%% extract TTL epochs/time stamps if recording contains any
if contains(nameGen{a,6},TTL)

    %storenameTTL = {'Epo1' 'PC0/'};
    storenameTTL = {'Epo1'};
    for k = 1:numel(storenameTTL)
      storename_TTL = storenameTTL{k};
      S4{k} = tdt2mat(tankdir, tankname, blockname, storename_TTL);
    end
end

 %%   
% store tdt2mat extracts in independent variable
if numel(storenames) == 1
    SynapseSignal = S{1};
end

if numel(storenames2) == 1
    SynapseReference1 = S2{1};
end

if numel(storenames3) == 1
    SynapseReference2 = S3{1};
end
%% Get raw traces for signal and reference channels

if strcmp(nameGen{a,10},FPNAPE)
    [raw470,raw565,ts,ts2] = extractRawTraces(S,S2);
elseif strcmp(nameGen{a,10},FPTwoMice1)
    [raw405,raw470,ts,ts2] = extractRawTraces(S,S2);
elseif strcmp(nameGen{a,10},FPTwoMice2)
    [raw405,raw470,ts,ts2] = extractRawTraces(S,S2);
elseif strcmp(nameGen{a,10},FPTwoMice1TEMP)
    [raw405,raw470,ts,ts2] = extractRawTraces(S,S2);
elseif strcmp(nameGen{a,10},FPTwoMice2TEMP)
    [raw405,raw470,ts,ts2] = extractRawTraces(S,S2);
else
    [raw470,raw405,raw565,ts,ts2,ts3] = extractRawTracesDUAL(S,S2,S3);
end

Dts = ts;
Fs = SynapseSignal.sampling_rate;

% adjust time stamp array from TTLs so that they start at correct time
if contains(nameGen{a,6},TTL)
    TTL_array1 = S4{1,1}.timestamps;
    ts_start = S{1,1}.timestamps(1);    %ts_start is the value where the recording started
    TTL_array = TTL_array1 - ts_start;  %subtract ts_start from all to get TTL event time relative to recording start
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subtracted 405 signal from 470 signal % raw final signal
if length(raw470) ~= length(raw405)
    if length(raw470) >= length(raw405)
        raw470 = raw470(1:length(raw405));
    elseif length(raw470) <= length(raw405)
        raw405 = raw405(1:length(raw470));
    end
end

if strcmp(nameGen{a,10},FPNAPE)
    rawsub565 = raw470-raw565;
elseif strcmp(nameGen{a,10},FPTwoMice1)
    rawsub405 = raw470-raw405;
elseif strcmp(nameGen{a,10},FPTwoMice2)
    rawsub405 = raw470-raw405;
elseif strcmp(nameGen{a,10},FPTwoMice1TEMP)
    rawsub405 = raw470-raw405;
elseif strcmp(nameGen{a,10},FPTwoMice2TEMP)
    rawsub405 = raw470-raw405;
else
    rawsub405 = raw470-raw405;
    rawsub565 = raw470-raw565;
end

if strcmp(nameGen{a,10},FPTwoMice1)
    dualcolor = 0;
elseif strcmp(nameGen{a,10},FPTwoMice2)
    dualcolor = 0;
elseif strcmp(nameGen{a,10},FPTwoMice1TEMP)
    dualcolor = 0;
elseif strcmp(nameGen{a,10},FPTwoMice2TEMP)
    dualcolor = 0;
else
    dualcolor = 1;
end

    % plot raw 405 ch, raw 470 ch and subtracted signal
    figure(1)
    plot(Dts,raw470,'b');
    hold on
    if dualcolor == 1
    plot(Dts,raw565,'r');
    plot(Dts,rawsub565,'m');
    end
    title('all channels (raw)')
    if ~strcmp(nameGen{a,10},FPNAPE)
        plot(Dts,raw405,'k');
        plot(Dts,rawsub405,'g');
    end
    
    if contains(nameGen{a,6},TTL)
        for k = 1:length(TTL_array)
            xline(TTL_array(k));
        end
    end
        
        
    xlabel('time(s)')
    ylabel('amplitude')
    xlim([1 ts(end)])
    %ylim([-30 210])
    hold off

    FigLocation = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/fiberphotometry_DATA/figures_rawtrace';    
    fullFileFig = fullfile(FigLocation,figname);
    saveas(gcf,fullFileFig,'epsc');
   






%% Save variables into a struct
if strcmp(nameGen{a,10},FPNAPE)
    FPvarnames = {'Dts' 'raw470' 'raw565' 'rawsub565'};
    FPvariables = [Dts raw470 raw565 rawsub565];
elseif strcmp(nameGen{a,10},FPTwoMice1)
    FPvarnames = {'Dts' 'raw405' 'raw470' 'rawsub405'};
    FPvariables = [Dts raw405 raw470 rawsub405];
elseif strcmp(nameGen{a,10},FPTwoMice2)
    FPvarnames = {'Dts' 'raw405' 'raw470' 'rawsub405'};
    FPvariables = [Dts raw405 raw470 rawsub405];
elseif strcmp(nameGen{a,10},FPTwoMice1TEMP)
    FPvarnames = {'Dts' 'raw405' 'raw470' 'rawsub405'};
    FPvariables = [Dts raw405 raw470 rawsub405];
elseif strcmp(nameGen{a,10},FPTwoMice2TEMP)
    FPvarnames = {'Dts' 'raw405' 'raw470' 'rawsub405'};
    FPvariables = [Dts raw405 raw470 rawsub405];
% elseif strcmp(nameGen{a,10},BruchasCartDual490)
%     FPvarnames = {'Dts' 'raw470' 'raw405' 'raw565' 'rawsub405'};
%     FPvariables = [Dts raw470 raw405 raw565 rawsub405];
else  
    FPvarnames = {'Dts' 'raw470' 'raw405' 'raw565' 'rawsub405' 'rawsub565'};
    FPvariables = [Dts raw470 raw405 raw565 rawsub405 rawsub565];
end

for i = 1:numel(FPvarnames)
    FPvalues.(string(FPvarnames(i))) = FPvariables(:,i);
end


%% Save time stamps as .txt file if TTLs were used for epochs
if contains(nameGen{a,6},TTL)
    %tsname_to_save = {'TTL_array'};
    %tsval_to_save = {TTL_array};
    timestamp_location = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/fiberphotometry_DATA/timestamps';
    TTL_location = append(timestamp_location,'/',TTLname);
    writematrix(TTL_array,TTL_location)
end

%% Save file as .mat file with specified filename

DataLocation = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/fiberphotometry_DATA/matfiles_extracted';    
fullFileData = fullfile(DataLocation,filename);
save(fullFileData,'Fs','FPvalues','FPvarnames');


%save(filename,'Fs','FPvalues','FPvarnames');

close all
clearvars -except a

end


%% Alternative DF/F methods

% GEf1 = polyfit(Dts,raw470,4);
% GEfitcurve1 = polyval(GEf1,Dts);
% GE470fit = raw470 - GEfitcurve1;
% adjusted470 = GE470fit;
% 
% GEf2 = polyfit(Dts,raw405,4);
% GEfitcurve2 = polyval(GEf2,Dts);
% GE405fit = raw405 - GEfitcurve2;
% 
% 
% % bls = polyfit(GE405fit,adjusted470,1);
% % Y_fit_all = bls(1).*GE405fit+bls(2);
% % Y_dF_all = adjusted470 - Y_fit_all;
% % dFF = (Y_dF_all)./(Y_fit_all);






