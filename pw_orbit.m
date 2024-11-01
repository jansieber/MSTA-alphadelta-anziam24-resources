function [iv,res]=pw_orbit(~,iv,u)
[     y01,      y02,      y03,      bstar,      jnew,      G,      alpha2,      alpha4]=deal(...
 u(iv.y01),u(iv.y02),u(iv.y03),u(iv.bstar),u(iv.jnew),u(iv.G),u(iv.alpha2),u(iv.alpha4));
[      ts1off,      ts1on,      ts2off,      ts2on,      T]=deal(...
  u(iv.ts1off),u(iv.ts1on),u(iv.ts2off),u(iv.ts2on),u(iv.T));
[      t1min,        y1_min,      y1dot_min]=deal(...
  u(iv.t1min),u(iv.y1_min),u(iv.y1dot_min));
% more formulas here
b=[1;bstar];
u=[jnew/G;jnew*alpha4/bstar];
y1per=yper(ts1on-ts1off,T,b(1),u(1));
y2per=yper(ts2on-ts2off,T,b(2),u(2));
y1_ts2off=pw_advance(y1per,     ts2off-ts1off,b(1),   0);
%y1_ts2on =pw_advance(y1per,     ts2on-ts1off, b(1),u(1));
y1_ts2on =pw_advance(y1per,     ts2on-(T+ts1off), b(1),u(1));

%y2_ts1off=pw_advance(y2per,   T+ts1off-ts2off,b(2),u(2));
y2_ts1off=pw_advance(y2per,   ts1off-ts2off,b(2),u(2));

%y2_ts1on =pw_advance(y2per,   T+ts1on-ts2off, b(2),   0);
y2_ts1on =pw_advance(y2per,   ts1on-ts2off, b(2),   0);

y1_ts1on =pw_advance(y1per,     ts1on-ts1off, b(1),   0);
%y1_ts1on =pw_advance(y1per,     ts1on-(T+ts1off), b(1),   0);


y1t1min  =pw_advance(y1_ts1on,  t1min-ts1on,  b(1),u(1));
y3=alpha2*jnew/G;
res=[
    y1_ts2off(1)-y02;...
    y1_ts2on(1)- y02;...
    y2_ts1off(1)-y3+y01;...
    y2_ts1on(1)- y3+y01;...
    y1t1min-[y1_min;y1dot_min]];
end
%%
function y=yper(toff2on,T,b,u)
M=@(t)pw_lin(t,b);
Id=eye(2);
y=(Id-M(T))\(Id-M(T-toff2on))*[u;0];
end
