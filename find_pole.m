function [p,i] = find_pole(A)
n = length(A);
p = [];
i = [];
for m=1:n-2
  if((A(m+1)>A(m))&&(A(m+1)>A(m+2)))
      i = [i,m+1];
      p = [p,A(m+1)];
  elseif((A(m+1)<A(m))&&(A(m+1)<A(m+2)))
      i = [i,m+1];
      p = [p,A(m+1)];
  end
end
end