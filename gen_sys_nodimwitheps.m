%% generate right-hand side for non dim Rit-Jansen model
%% load dim param 
 dim=load('C:\Users\hm672\OneDrive - University of Exeter\Documents\MATLAB\gen-nondimritjansen Coco functions\dimeparam.mat');
parnames={'jnew','G','P','alpha1','alpha2','alpha3','alpha4','b_star', ...
    'c','beta','epsilon',...
    'x01','x02','x03','a1','a2','a3'};
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either j or
% par(1) to refer to parameter j
syms(parnames{:});       % create symbols for G,...
par=cell2sym(parnames);  % now j is par(1) etc
y=sym('y',[3,1]);        % create symbols for dependent variables y
yp=sym('yp',[3,1]);      % and derivatives
syms z;                   % add new var for y3-y2
%% useful functions and constants
% Dimensionless func. 
newsigm=@(x,x0,epsilon)1./(1+exp(((x0)-x)./(epsilon)));

%% r.h.s
% jnew=2
dt(1:3)=yp;
dtt(1)= (jnew/G)*newsigm(y(3)-y(2),x01,epsilon/a1)        -2*yp(1)        -y(1);
dtt(2)= jnew*(b_star*alpha4*newsigm(y(1),x02,epsilon/a2)) -2*b_star*yp(2)-b_star^2*y(2);
dtt(3)=(P/G)+(jnew/G)*(alpha2*newsigm(y(1),x03,(epsilon/a3))) -2*yp(3)       -y(3);
dz=dt(3)-dt(2)-(z-(y(3)-y(2)));
f=[dt(:);dtt(:);dz];
%% generate r.h.s. file and derivative
sco_sym2funcs(f,...           % symbolic expression for f
    {[y;yp;z],par},...          % which symbols are in which inputs of f
    {'x','p'},...             % names for inputs of f
    'vector',[1,1],...        % are inputs scalar or vectors (default: vectors)
    'filename','symnondritjan_eps');% filename for result
