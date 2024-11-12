clear;
varnames={'bstar','jnew','G','alpha2','alpha4',...
    'a1','a2','a3','y01','y02','y03','y1','epsilon'};%% Create symbols for parameters
syms(varnames{:});       % create symbols for G,bstar,...
u=cell2sym(varnames);
cl=[varnames;num2cell(1:length(varnames))];
iv=struct(cl{:}); 
%% useful functions and constants
S = @(x, x0, epsilon) 1 ./ (1 + exp((x0 - x) ./ (epsilon)));

% Define Y2 and Y3 functions
y2 =  (alpha4 * jnew / bstar) * S(y1, y02, epsilon / a2);
y3 =  (alpha2 * jnew / G) * S(y1, y03, epsilon / a3); 
y = S((y3 - y2),y01 , epsilon / a1);
ydiff=simplify(diff(y,y1))
output=matlabFunction([y,ydiff],'file','functionofu','Vars',{u})