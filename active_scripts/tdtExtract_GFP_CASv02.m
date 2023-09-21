% Code written by Carrie Stine, parts adapted from code written by Christian Pedersen
% Michael Bruchas Lab - UW

% Extracts photometry data from TDT Data Tank and outputs time series as
% .mat file

% Updated 9/21/23 to V2, updated file locations to new GitHub repository

% REQUIRES tdt2mat.m file to be in same folder (filepath)

%% Reset MatLab workspace - clears all variables and the command window

clear variables;  % clear all variables
close all;  % close all open graphs
clc   % clear command window


%% 
A = [126];
for a1 = A %change this range to which rows of your nameGenerator file you want to run
   
 
    %%%%%%%%%%%%%%%%%%%%%%  EDIT THESE FIELDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
% point to tanks (enter file path)
path_to_data = '/Users/castine/Dropbox (Bruchas Lab)/Carrie Data/Bruchas Lab/pnVTA nociceptin/DATA/NOPLight validation/Tanks'; %change to point to where your data is

% set up Tank name, variables to extract
tankdir = (path_to_data);

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
%Bruchas or NAPE cart:
%storenames = {'470A'}; % name of stores to extract from TDT (usu. 4-letter code) 
%LMag is the demodulated data, may also have other timestamps etc

%storenames2 = {'405A'};

FPstaticA = 'StaticA';
FPstaticB = 'StaticB';
FPstaticC = 'StaticC';
FPstaticD = 'StaticD';
FPstaticC_TEMP = 'StaticC_TEMP';
FPstaticD_TEMP = 'StaticD_TEMP';
MixedLED490 = 'BruchasCart_MixedLED490';
BruchasCart490 = 'BruchasCart490';
BruchasCart470B = 'BruchasCart470B';
BruchasCart435470 = 'BruchasCart435470';
TTL = 'TTL';

if strcmp(nameGen{a,10},FPstaticA)
    storenames = {'470A'};
    storenames2 = {'405A'};
elseif strcmp(nameGen{a,10},FPstaticB)
    storenames = {'470B'};
    storenames2 = {'405B'};
elseif strcmp(nameGen{a,10},FPstaticC)
    storenames = {'465A'};
    storenames2 = {'405A'};
elseif strcmp(nameGen{a,10},FPstaticD)
    storenames = {'465C'};
    storenames2 = {'405C'};   
elseif strcmp(nameGen{a,10},FPstaticC_TEMP)
    storenames = {'405C'};
    storenames2 = {'465C'};   
elseif strcmp(nameGen{a,10},FPstaticD_TEMP)
    storenames = {'405B'};
    storenames2 = {'465B'};     
elseif strcmp(nameGen{a,10},MixedLED490)
    storenames = {'490B'};
    storenames2 = {'435B'};
elseif strcmp(nameGen{a,10},BruchasCart490)
    storenames = {'490A'};
    storenames2 = {'435A'};
elseif strcmp(nameGen{a,10},BruchasCart470B)
    storenames = {'470B'};
    storenames2 = {'405B'};
elseif strcmp(nameGen{a,10},BruchasCart435470)
    storenames = {'470A'};
    storenames2 = {'435A'};
else
    storenames = {'470A'};
    storenames2 = {'405A'};
end

% % Static system Box C:
% storenames = {'465C'}; % name of stores to extract from TDT (usu. 4-letter code) 
% %LMag is the demodulated data, may also have other timestamps etc
% 
% storenames2 = {'405C'};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract
for k = 1:numel(storenames)
  storename = storenames{k};
  S{k} = tdt2mat(tankdir, tankname, blockname, storename); %#ok<SAGROW>
end

for k = 1:numel(storenames2)
  storename2 = storenames2{k};
  S2{k} = tdt2mat(tankdir, tankname, blockname, storename2); %#ok<SAGROW>
end

if numel(storenames) == 1
    SynapseSignal = S{1};
end

if numel(storenames2) == 1
    SynapseReference = S2{1};
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

%% Get raw traces for signal and reference channels

[raw470,raw405,ts,ts2] = extractRawTraces(S,S2);
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
rawsub405 = raw470-raw405;



    % plot raw 405 ch, raw 470 ch and subtracted signal
    figure(1)
    plot(Dts,raw470,'b');
    hold on
    plot(Dts,raw405,'k');
    title('all channels (raw)')
    plot(Dts,rawsub405,'g');
    
    if contains(nameGen{a,6},TTL)
        for k = 1:length(TTL_array)
            xline(TTL_array(k));
        end
    end
    
    xlabel('time(s)')
    ylabel('amplitude')
    xlim([1 ts(end)])
    %ylim([-30 50])
    hold off

FigLocation = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository/fiberphotometry_DATA/figures_rawtrace';    
fullFileFig = fullfile(FigLocation,figname);
saveas(gcf,fullFileFig,'epsc');





%% Save variables into a struct
FPvarnames = {'Dts' 'raw470' 'raw405' 'rawsub405'};

FPvariables = [Dts raw470 raw405 rawsub405];

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

