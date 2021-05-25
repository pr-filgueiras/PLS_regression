function [y_hat,D]= savgol(y,width,order,deriv)
%SAVGOL Savitsky-Golay smoothing and differentiation.
%  Inputs are the matrix of ROW vectors to be smoothed (y),
%  and the optional variables specifying the number of points in
%  filter (width), the order of the polynomial (order), and the
%  derivative (deriv). The output is the matrix of smoothed
%  and differentiated ROW vectors (y_hat) and the matrix of 
%  coefficients (cm) which can be used to create a new smoothed/
%  differentiated matrix, i.e. y_hat = y*cm. If number of points,
%  polynomial order and derivative are not specified,
%  they are set to 15, 2 and 0, respectively.
%
%  Example: if y is a 5 by 100 matrix then savgol(y,11,3,1)
%  gives the 5 by 100 matrix of first-derivative row vectors
%  resulting from a 11-point cubic Savitzky-Golay smooth of
%  each row of y.
%
%I/O format is: [y_hat,cm] = savgol(y,width,order,deriv);
%
%See also: MSCORR, SAVGOLCV, SGDEMO, STDFIR, BASELINE, LAMSEL


% Sijmen de Jong Unilever Research Laboratorium Vlaardingen Feb 1993
% Modified by Barry M. Wise 5/94
%         ***   Further modified, 1998-03, Martin Andersson
%         ***   Adjusting the calcn. of the bulk data.
%         ***   Based on calcn. of a sparse derivative matrix (D)

[m,n] = size(y);
y_hat = y;
% set default values: 15-point quadratic smooth
if nargin<4
  deriv= 0;
  disp('  '), disp('Derivative set to zero') 
end
if nargin<3
  order= 2; 
  disp('  '), disp('Polynomial order set to 2')
end
if nargin<2
  width=min(15,floor(n/2)); 
  s = sprintf('Width set to %g',width);
  disp('  '), disp(s)  
end
% In case of input error(s) set to reasonable values
w = max( 3, 1+2*round((width-1)/2) );
if w ~= width
  s = sprintf('Width changed to %g',w);
  disp('  '), disp('Width must be >= 3 and odd'), disp(s)
end
o = min([max(0,round(order)),5,w-1]);
if o ~= order
  s = sprintf('Order changed to %g',o); disp('  ')
  disp('Order must be <= width -1 and <= 5'), disp(s)
end
d = min(max(0,round(deriv)),o);
if d ~= deriv
  s = sprintf('Derivative changed to %g',d); disp('  ')
  disp('Deriviative must be <= order'), disp(s)
end
p = (w-1)/2;
% Calculate design matrix and pseudo inverse
x = ((-p:p)'*ones(1,1+o)).^(ones(size(1:w))'*(0:o));
weights = x\eye(w);
% Smoothing and derivative for bulk of the data
coeff=prod(ones(d,1)*[1:o+1-d]+[0:d-1]'*ones(1,o+1-d,1),1);
D=spdiags(ones(n,1)*weights(d+1,:)*coeff(1),p:-1:-p,n,n);
% Smoothing and derivative for tails 
w1=diag(coeff)*weights(d+1:o+1,:);
D(1:w,1:p+1)=[x(1:p+1,1:1+o-d)*w1]'; 
D(n-w+1:n,n-p:n)=[x(p+1:w,1:1+o-d)*w1]';
% Operate on y using the filtering/derivative matrix, D
y_hat=y*D;
