function y=pw_advance(y0,tc,b,u)
for i=length(tc):-1:1
    M=pw_lin(tc(i),b);
    y(:,i)=M*y0+(eye(2)-M)*[u;0];
end
end
