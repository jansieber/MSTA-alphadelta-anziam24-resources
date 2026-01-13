%% generate right-hand side for FN-model
clear
parnames={'a','epsilon'};
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either A or
% par(1) to refer to parameter A
syms(parnames{:});       % create symbols for A,B,...
par=cell2sym(parnames);  % now A is par(1) etc
y=sym('y',[2,1]);        % create symbols for dependent variables y
%% useful functions and constants
nul=@(x)(x.^3/3)-x; 
%% r.h.s  % x',y'
yp(1)= (y(2) -nul(y(1)))./epsilon;
yp(2)=a-y(1);
f=[yp(1);yp(2)];
%% generate r.h.s. file and derivative
sco_sym2funcs(f,...           % symbolic expression for f
    {y,par},...          % which symbols are in which inputs of f
    {'x','p'},...             % names for inputs of f
    'vector',[1,1],...        % are inputs scalar or vectors (default: vectors)
    'filename','sym_FN');% filename for result
%% run this file and then 
% sym_FN Automatically generated with matlabFunction