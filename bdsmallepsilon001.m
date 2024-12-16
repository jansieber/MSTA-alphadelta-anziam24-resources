% run startup COCO function
% run file gen_sys_nodimwitheps

%% Load path
% clear
cocofolder=fileparts(which('coco'));
symcocopath=[cocofolder,'/../../contributed/symcoco/toolbox'];
addpath(symcocopath); % path of symcoco routines
parnames={'jnew','G','P','alpha1','alpha2','alpha3','alpha4','b_star','c','beta',...
    'epsilon','x01','x02','x03','a1','a2','a3'};
format compact
cl=[parnames;num2cell(1:length(parnames))];
ip=struct(cl{:});                % automaticaly number parameters
F=sco_gen(@symnondritjan_eps);       % F and F2 are the same
funcs = {F(''),F('x'),F('p')};   % r.h.s, derivative rt x and derivative wrt p
%% initial guess for parameters
dimpar=load([pwd(),filesep,'dimeparam.mat']);% dowonload param scribt file 
beta=2*dimpar.a/(dimpar.B*dimpar.r*dimpar.C*2*dimpar.e0);%0.024; 
alpha1=1;alpha2=0.8; alpha3=0.25;
jnew=2;
G=0.5;
c=dimpar.r*dimpar.v0;
x01=c*beta;
x02=c*beta/alpha3;                                                                                                                                                          
x03=c*beta/alpha1;
a1=1;
a2=alpha3;
a3=alpha1;
par([ip.x01;ip.x02;ip.x03;ip.a1;ip.a2;ip.a3;ip.beta;ip.epsilon])=...
    [x01;     x02;    x03; a1;   a2;    a3;   beta;  0.001];
par([ip.jnew;ip.G; ip.P;ip.alpha1;ip.alpha2;ip.alpha3;ip.alpha4;ip.b_star;ip.c ])=...
   [jnew;     G ;    0;        alpha1; 0.8;  alpha3;    0.25 ; 0.5; c];
dim=7;
x0=[1.14279;0.815087;0.9425519;0.0265967;0.068872;6.1237e-05;0.153799];%bifurc
% if we start from zero ini condition 
% x0=zeros(dim,1); % then there no bifurc
% only zero xeq 
prob = coco_prob();
prob = ode_isol2ep(prob, '', funcs{:}, x0, parnames, par);
prob = coco_set(prob, 'cont', 'PtMX',  [50  60000],'h0',1e-4,'h_min',1e-6,'h_max',0.05,'NPR',500);
fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
  'ep_run');
coco(prob, 'smalleps=001', [], 1, 'G', [0.3 20]);


%% max and min
bd_ep_P_0=coco_bd_read('smalleps=001');
HBlab = coco_bd_labs(bd_ep_P_0, 'HB');
prob = coco_prob();
prob = coco_set(prob, 'coll', 'MXCL',false); % remove discretization error limit
prob = ode_HB2po(prob, '', 'smalleps=001', HBlab(1));
prob = coco_set(prob, 'cont', 'PtMX',[200,0], 'NAdapt', 5,...
   'h_max',1e5,'NPR',100,'MaxRes',100,'norm',inf);
coco(prob, 'po_smalleps001UZ', [], 1, 'G', [0.5 4]);
%%
% data for canrad near Hopf bifurcation G1 

% define user points near Hopf bifurcation G1 
% prob = coco_add_event(prob, 'uz', 'G', 1.4821:0.0001:1.5);
% prob = coco_set(prob, 'cont', 'PtMX',[200,0], 'NAdapt', 5,...
%    'h_max',1e5,'NPR',100,'MaxRes',100,'norm',inf);
% coco(prob, 'po_Canardeps001UZ', [], 1, 'G', [0.5 2]);
%%
