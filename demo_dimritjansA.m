%% Load path (run by coco March 2023)
%% run file gen_sym_ritjansen
%%
clear
  cocofolder=fileparts(which('coco'));
  symcocopath=[cocofolder,'/../../contributed/symcoco/toolbox'];
  addpath(symcocopath); % path of symcoco routines
 parnames={'A','B','a','b','C1','e0','r','v0','P'};
format compact
cl=[parnames;num2cell(1:length(parnames))];
ip=struct(cl{:});                % automaticaly number parameters
F=sco_gen(@sym_ritjansen);       % F and F2 are the same
funcs = {F(''),F('x'),F('p')};   % r.h.s, derivative rt x and derivative wrt p
%% bd in paramter A
% initial guess for parameters
par([ip.A;ip.B;ip.a;ip.b;ip.C1;ip.e0; ip.r;ip.v0;ip.P])=...
    [3.25;  22; 100;  50;  135;  2.5;  0.56;    6;   0]; 
dim=7;
prob = coco_prob();
prob = ode_isol2ep(prob, '', funcs{:}, zeros(dim,1), parnames, par);
prob = coco_set(prob, 'cont', 'PtMX', 5000,'h_max',0.1,'NPR',100,'MaxRes',100);
fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
  'ep_run');
% filename= 'ep_inputp0'
bd_ep_P  = coco(prob, 'ep_inputp0', [], 1, 'A', [0.1 20]);

%% branch off at Hopf
bd_ep_P_0=coco_bd_read('ep_inputp0');
HBlab = coco_bd_labs(bd_ep_P_0, 'HB');
prob = coco_prob();
prob = coco_set(prob, 'coll', 'MXCL',false); % remove discretization error limit
prob = ode_HB2po(prob, '', 'ep_inputp0', HBlab(3));
prob = coco_set(prob, 'cont', 'PtMX',[700,0], 'NAdapt', 5,...
   'h_max',1e3,'NPR',10,'MaxRes',100,'norm',inf);
coco(prob, 'po_inputp0', [], 1, 'A', [0.5 20]);
%%
%p=120

par([ip.A;ip.B;ip.a;ip.b;ip.C1;ip.e0; ip.r;ip.v0;ip.P])=...
    [0.5;  22; 100;  50;  135;  2.5;  0.56;    6;   120]; 
% A=3.5 output mx 
dim=7;
prob = coco_prob();
prob = ode_isol2ep(prob, '', funcs{:}, zeros(dim,1), parnames, par);
prob = coco_set(prob, 'cont', 'PtMX', 5000,'h_max',0.1,'NPR',100,'MaxRes',100);
fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
  'ep_run');
% filename= 'ep_inputp0'
bd_ep_P  = coco(prob, 'ep_inputp120', [], 1, 'A', [0.5 20]);
%%

bd_ep_P_0=coco_bd_read('ep_inputp120');
HBlab = coco_bd_labs(bd_ep_P_0, 'HB');
prob = coco_prob();
prob = coco_set(prob, 'coll', 'MXCL',false); % remove discretization error limit
prob = ode_HB2po(prob, '', 'ep_inputp120', HBlab(3));
prob = coco_set(prob, 'cont', 'PtMX',[700,0], 'NAdapt', 5,...
   'h_max',1e3,'NPR',10,'MaxRes',100,'norm',inf);
coco(prob, 'po_inputp120', [], 1, {'A','po.period'}, [0.5 20]);
% this data used to produce time profile 



