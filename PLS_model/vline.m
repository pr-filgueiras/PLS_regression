function h = vline(x,lc)
[m,n] = size(x);
if m>1&n>1
  error('Error - input must be a scaler or vector')
elseif n>1
  x   = x';
  m   = n;
end

v     = axis;
if ishold
  for ii=1:m
    h(ii) = plot([1 1]*x(ii,1),v(3:4),lc);
  end
else
  hold on
  for ii=1:m
    h(ii) = plot([1 1]*x(ii,1),v(3:4),lc);
  end
  hold off
end

if nargout == 0;
  clear h
end