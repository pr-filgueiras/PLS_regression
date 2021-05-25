function h = hline(y,lc)
[m,n] = size(y);
v = axis;
if ishold
  for ii=1:m
    h(ii) = plot(v(1:2),[1 1]*y(ii,1),lc);
  end
else
  hold on
  for ii=1:m
    h(ii) = plot(v(1:2),[1 1]*y(ii,1),lc);
  end
  hold off
end

if nargout == 0;
  clear h
end
