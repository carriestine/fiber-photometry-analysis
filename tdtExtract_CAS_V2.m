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

%%%%%%%%%%%%%%%%%%%%%%  EDIT THESE FIELDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
% point to tanks (enter file path)
path_to_data = '/Users/castine/Dropbox (Bruchas Lab)/TEMPORARY PHOTOMETRY/Carrie_analysis_repository'; %change to point to where your data is

% set up Tank name, variables to extract
tankdir = [path_to_data];


nameGenerator = readtable('nameGeneratorVTAGCAMP.csv'); %enter name of your excel sheet containing naming info
nameGen = table2cell(nameGenerator);

plotFigures = 0; %change to 1 if you want to see the plots
%%
for a = 1:5 %change this range to which rows of your nameGenerator file you want to run
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    tankname = nameGen{a,1};
    
    blockname = nameGen{a,2};
    
    tempID = string(nameGen(a,:));
    
    filename1 = strjoin(tempID(3:5),'');
    filename2 = append(filename1,'.',"mat");
    filename = char(filename2);





%%
%Bruchas cart:
storenames = {'470A'}; % name of stores to extract from TDT (usu. 4-letter code) 
%LMag is the demodulated data, may also have other timestamps etc

storenames2 = {'405A'};


% % Static system Box C:
% storenames = {'465C'}; % name of stores to extract from TDT (usu. 4-letter code) 
% %LMag is the demodulated data, may also have other timestamps etc
% 
% storenames2 = {'405C'};



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

if numel(storenames) == 1
    SynapseSignal = S{1};
end

if numel(storenames2) == 1
    SynapseReference = S2{1};
end

%% Get raw traces for signal and reference channels

[raw470,raw405,ts,ts2] = extractRawTraces(S,S2);
Dts = ts;
Fs = SynapseSignal.sampling_rate;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subtracted 405 signal from 470 signal % raw final signal
rawsub405 = raw470-raw405;


if plotFigures == 1
    % plot raw 405 ch, raw 470 ch and subtracted signal
    figure(1)
    plot(Dts,raw470,'b');
    hold on
    plot(Dts,raw405,'r');
    title('both signals')
    plot(Dts,rawsub405,'g');
    xlabel('time(s)')
    ylabel('amplitude')
    xlim([1 ts(end)])
    ylim([-30 210])
    hold off

end


%% fit curve to adjust for signal decay

[adjusted470, adjusted470sub405,fitcurve470,fitcurve405,GE405fit2] = adjustDecay(raw470, raw405, Dts);



%% Calculate dF/F

[dFF470] = getdFF(adjusted470);

%dFF470sub405 = (adjusted470-GE405fit2);

[dFF470sub405] = getdFF(adjusted470,GE405fit2);


%% Calculate Z Score
[Z470] = getZscore(adjusted470);
[ZdFF470] = getZscore(dFF470);

[Z470sub405] = getZscore(adjusted470sub405);
[ZdFF470sub405] = getZscore(dFF470sub405);

%%

if plotFigures == 1
    % plot baseline corrected reference
    figure(101);
    plot(Dts,raw405,Dts,fitcurve405,Dts,GE405fit2);
    title('Signal with corrected baseline')
    xlabel('time(s)')
    ylabel('raw F')
    
    % plot baseline corrected signal
    figure(102);
    plot(Dts,raw470,Dts,fitcurve470,Dts,adjusted470);
    title('Signal with corrected baseline')
    xlabel('time(s)')
    ylabel('raw F')
    
    % plot baseline corrected (signal - reference) vs raw (signal - reference)
    figure(103);
    plot(Dts,rawsub405,Dts,adjusted470sub405);
    title('Signal - reference with corrected baseline')
    xlabel('time(s)')
    ylabel('raw F')
    
    %%
    % plot DFF baseline-corrected signal
    figure(201);
    plot(Dts,raw470,'k')
    hold on
    plot(Dts,adjusted470,'b')
    plot(Dts,dFF470,'g')
    title('DF/F signal')
    xlabel('time(s)')
    ylabel('∆F/F')
    hold off
    
    % plot DFF baseline-corrected signal
    figure(202);
    plot(Dts,rawsub405,'k')
    hold on
    plot(Dts,adjusted470sub405,'b')
    plot(Dts,dFF470sub405,'g')
    title('DF/F (signal - reference)')
    xlabel('time(s)')
    ylabel('∆F/F')
    hold off
    
    %%
    
    % plot Z-scored DFF baseline-corrected signal
    figure(301);
    plot(Dts,Z470,'g');
    hold on
    plot(Dts,ZdFF470,'r');
    title('Z-scored signal')
    xlabel('time(s)')
    ylabel('Z score')
    hold off
    
    figure(302)
    plot(Dts,Z470sub405,'g');
    hold on
    plot(Dts,ZdFF470sub405,'r');
    title('Z-scored(signal - reference)')
    xlabel('time(s)')
    ylabel('Z score')
    hold off
end


%% Save variables into a struct
FPvarnames = {'Dts' 'raw470' 'raw405' 'rawsub405' 'adjusted470' 'adjusted470sub405' 'dFF470' 'dFF470sub405' 'Z470' 'Z470sub405' 'ZdFF470' 'ZdFF470sub405'};

FPvariables = [Dts raw470 raw405 rawsub405 adjusted470 adjusted470sub405 dFF470 dFF470sub405 Z470 Z470sub405 ZdFF470 ZdFF470sub405];

for i = 1:numel(FPvarnames)
    FPvalues.(string(FPvarnames(i))) = FPvariables(:,i);
end

%% Save file as .mat file with specified filename

save(filename,'Fs','FPvalues');

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






