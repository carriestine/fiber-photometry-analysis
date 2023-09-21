
function [signal,reference1,reference2,ts,ts2,ts3] = extractRawTracesDUAL(S,S2,S3)
% Massage data and get time stamps

LMag = S{1}; %add more if you extracted more stores above
% LMag2 = S{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani1 = LMag.channels==1;
% chani2 = LMag.channels==2;

LMag2 = S2{1}; %add more if you extracted more stores above
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani21 = LMag2.channels==1;
% chani22 = LMag2.channels==2;

LMag3 = S3{1}; %add more if you extracted more stores above
% LMag2 = S{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani31 = LMag3.channels==1;
%chani32 = LMag3.channels==2;

% Get LMag data as a vector (repeat for each channel)
rawdat1 = LMag.data(chani1,:);
signal = reshape(rawdat1', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all channels)
ts = LMag.timestamps(chani1);
%t_rec_start = ts(1);

ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
ts = bsxfun(@plus, ts(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
ts = reshape(ts',[],1);

%%%%%%%%%%%%%%%%%%


dat2 = LMag2.data(chani21,:);
reference1 = reshape(dat2', [],1); % unwrap data from m x 256 array


% Get LMag timestamps (use chani1 - timestamps should be the same for all channels)
ts2 = LMag2.timestamps(chani21);
%t_rec_start2 = ts2(1);

ts2 = ts2-ts2(1); % convert from Unix time to 'seconds from block start'
ts2 = bsxfun(@plus, ts2(:), (0:LMag2.npoints-1)*(1./LMag2.sampling_rate));
ts2 = reshape(ts2',[],1);



dat3 = LMag3.data(chani31,:);
reference2 = reshape(dat3', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts3 = LMag3.timestamps(chani31);
%t_rec_start3 = ts3(1);

ts3 = ts3-ts3(1); % convert from Unix time to 'seconds from block start'
ts3 = bsxfun(@plus, ts3(:), (0:LMag3.npoints-1)*(1./LMag3.sampling_rate));
ts3 = reshape(ts3',[],1);








