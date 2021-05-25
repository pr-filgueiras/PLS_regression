function [RMSE,R2,bias] = rmse(y1,y2,nlv)
% Routine for computing merit figures as ASTM E 1655-05
% RMSE : root mean square error;
% R2   : determination coefficient; 
% bias : bias.
%
% input:
%  y1  : reference values;
%  y2  : predict values;
% nlv  : number of PLS latent variables for calibration error (RMSEC) 
%
% Paulo R. Filgueiras - 30/09/2013
%

if length(y1)~=length(y2)
    disp('')
    disp('size vectors differents')
    disp('')
return
end

if nargin==2 
    RMSE=sqrt((sum((y1-y2).^2))./length(y1));
    R2 = corrcoef(y1,y2); R2=R2(1,2)^2;
    bias = sum(y1-y2)/length(y1);
else
    RMSE=sqrt((sum((y1-y2).^2))./(length(y1)-nlv-1));
    R2 = corrcoef(y1,y2); R2=R2(1,2)^2;
    bias = sum(y1-y2)/length(y1);
end
