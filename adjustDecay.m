function [adjusted470,adjusted470sub405,GEfitcurve1,GEfitcurve2,GE405fit] = adjustDecay(raw470, raw405, Dts)

GEf1 = polyfit(Dts,raw470,4);
GEfitcurve1 = polyval(GEf1,Dts);
GE470fit = raw470 - GEfitcurve1;
adjusted470 = GE470fit;

GEf2 = polyfit(Dts,raw405,4);
GEfitcurve2 = polyval(GEf2,Dts);
GE405fit = raw405 - GEfitcurve2;

adjusted470sub405 = GE470fit - GE405fit;

% Robust non negative linear regression
    fitdata = fit(GE405fit,GE470fit,fittype('poly1'),'Robust','on');
  % Align reference trace to signal using the fit
    GE405fit2 = fitdata(GE405fit);
    
adjusted470sub405 = GE470fit - GE405fit2;

