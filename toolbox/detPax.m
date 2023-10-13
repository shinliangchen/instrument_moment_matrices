function y = detPax(nx,na)

% table of deterministic probability distributions for single side

for lambda = 1:na^nx
    
    v(lambda,:) = dec2base(lambda-1,na,nx)-'0';
    
end

for lambda = 1:na^nx
for a=0:na-1
for x=0:nx-1

  y(lambda,a+1,x+1) = kronDel(a,v(lambda,x+1));
  %disp(y(lambda,a+1,x+1))

end
end
end

end