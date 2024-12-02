clear
parnames={'A','B','a','b','C1','e0','r','v0','P'};
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either A or
% par(1) to refer to parameter A
syms(parnames{:});       % create symbols for A,B,...
par=cell2sym(parnames);  % now A is par(1) etc
y=sym('y',[3,1]);        % create symbols for dependent variables y
yp=sym('yp',[3,1]);      % and derivatives
syms z;                  % add new var for y3-y2
%% useful functions and constants
sigm=@(x)2*e0/(1+exp(r.*(v0-x))); 
C2=0.8*C1;
C3=0.25*C1;
C4=0.25*C1;
%% r.h.s
dt(1:3)=yp;
dtt(1)=         A*a*sigm(y(3)-y(2))-2*a*yp(1)-a^2*y(1);
dtt(2)=      B*b*C4*sigm(C3*y(1))  -2*b*yp(2)-b^2*y(2);
dtt(3)=A*a*P+a*A*C2*sigm(C1*y(1))  -2*a*yp(3)-a^2*y(3);
dz=yp(3)-yp(2)-(z-(y(3)-y(2)));
f=[dt(:);dtt(:);dz];
%% generate r.h.s. file and derivative
sco_sym2funcs(f,...           % symbolic expression for f
    {[y;yp;z],par},...          % which symbols are in which inputs of f
    {'x','p'},...             % names for inputs of f
    'vector',[1,1],...        % are inputs scalar or vectors (default: vectors)
    'filename','sym_ritjansen');% filename for result
