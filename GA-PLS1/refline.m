function h = refline(slope,intercept)
%REFLINE Add a reference line to a plot.
%   REFLINE(SLOPE,INTERCEPT) adds a line with the given SLOPE and
%   INTERCEPT to the current figure.
%
%   REFLINE(SLOPE) where SLOPE is a two element vector adds the line
%        y = SLOPE(2) + SLOPE(1)*x 
%   to the figure. (See POLYFIT.)
%
%   H = REFLINE(SLOPE,INTERCEPT) returns the handle to the line object
%   in H.
%
%   REFLINE with no input arguments superimposes the least squares line on 
%   the plot based on points recognized by LSLINE.
%
%   See also POLYFIT, POLYVAL, LSLINE.   

%   Copyright 1993-2008 The MathWorks, Inc. 


if nargin == 0
   hh = lsline;
   if nargout >0
       h = hh;
   end
   return;
end

if nargin == 1
   if max(size(slope)) == 2
      intercept=slope(2);
      slope = slope(1);
   else
      intercept = 0;
   end
end

xlimits = get(gca,'Xlim');
ylimits = get(gca,'Ylim');

np = get(gcf,'NextPlot');
set(gcf,'NextPlot','add');

xdat = xlimits;
ydat = intercept + slope.*xdat;
maxy = max(ydat);
miny = min(ydat);

if maxy > ylimits(2)
  if miny < ylimits(1)
     set(gca,'YLim',[miny maxy]);
  else
     set(gca,'YLim',[ylimits(1) maxy]);
  end
else
  if miny < ylimits(1)
     set(gca,'YLim',[miny ylimits(2)]);
  end
end

if nargout == 1
   h = line(xdat,ydat);
   set(h,'LineStyle','-');
else
   hh = line(xdat,ydat);
   set(hh,'LineStyle','-');
end

set(gcf,'NextPlot',np);
