%% run file gen_sys_nodimwitheps
clear
parnames={'jnew','G','P','alpha1','alpha2','alpha3','alpha4','b_star','c','beta',...
    'epsilon','x01','x02','x03','a1','a2','a3'};
format compact
cl=[parnames;num2cell(1:length(parnames))];
ip=struct(cl{:});                % automaticaly number parameters
F=sco_gen(@symnondritjan_eps);       % F and F2 are the same
funcs = {F(''),F('x'),F('p')};   % r.h.s, derivative rt x and derivative wrt p
%% initial guess for parameters
dimpar=load([pwd(),filesep,'dimeparam.mat']);
beta=2*dimpar.a/(dimpar.B*dimpar.r*dimpar.C*2*dimpar.e0);% 0.024;
alpha1=1; alpha3=0.25;
jnew=2;
% initial value
G=0.5;
bstar=0.5;
c=dimpar.r*dimpar.v0;
% here call y01=x01
x01=c*beta;
x02=c*beta/alpha3;
x03=c*beta/alpha1;
a1=1;
a2=alpha3;

a3=alpha1;
par([ip.x01;ip.x02;ip.x03;ip.a1;ip.a2;ip.a3;ip.beta;ip.epsilon])=...
    [x01;     x02;    x03; a1;   a2;    a3;   beta;  0.024];
par([ip.jnew;ip.G; ip.P;ip.alpha1;ip.alpha2;ip.alpha3;ip.alpha4;ip.b_star;ip.c ])=...
    [jnew;     G ;    0;        alpha1; 0.8;  alpha3;    0.25 ;bstar; c];
dim=7;
% choose non zero xeq for epsilon=0.024
x0=[1.14279;0.815087;0.9425519;0.0265967;0.068872;6.1237e-05;0.153799];
%x0=zeros(dim,1);
prob = coco_prob();
prob = ode_isol2ep(prob, '', funcs{:}, x0, parnames, par);
prob = coco_set(prob, 'cont', 'PtMX', 5000,'h_max',0.01,'h_min',1e-6,'NPR',100);

fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
    'ep_run');

%epsilon=0.024

coco(prob, 'bdeps0244', [], 1, 'G', [0.5 20]);
%%
bdep=coco_bd_read('bdeps0244');
HBlab = coco_bd_labs(bdep, 'HB');
prob = coco_prob();
prob = ode_HB2po(prob, '', 'bdeps0244', HBlab(1));
fprintf(...
    '\n Run=''%s'': Continue periodic orbits from point %d in run ''%s''.\n', ...
    'po_run', HBlab(1), 'ep_run');

prob = coco_add_event(prob, 'UZ', 'po.period', 20:20:500);

prob = coco_set(prob, 'cont', 'PtMX',[2000,0], 'NAdapt', 2,...
    'h0',1e-1,'h_min',1e-5 ,'NPR',10,'MaxRes',100,'norm',inf);

bd_Ho2po_a2_alph3 = coco(prob, 'POinput0', [], 1, {'G','po.period'}, {[0.5 20] [0 510]});
%%
%track Hopf and SN
%continue family of Hopf bifuraction by varying G,bstar
%%
old_name='bdeps0244';
bd_ep=coco_bd_read(old_name);

HBlab = coco_bd_labs(bd_ep, 'HB');
SNlab = coco_bd_labs(bd_ep, 'SN');

i=1;
runid = sprintf('HB_Gbstar%d_run', i);

prob = coco_prob();
prob = ode_ep2HB(prob, '', old_name, HBlab(i));

prob = coco_add_event(prob, 'UZ', 'b_star', 0.2:0.001:0.5);%linspace(0.1,0.5,20)
[data, uidx] = coco_get_func_data(prob, 'ep.hb', 'data', 'uidx');
prob = coco_add_pars(prob, 'om_squared', uidx(data.ep_hb.k_idx), {'om2'});
fprintf(...
    '\n Run=''%s'': Continue family of HB from point %d in run ''%s''.\n', ...
    runid, i, 'HB_run');
coco(prob,runid , [], 1, {'G' 'b_star' 'om2'}, {[0.3 3.5] [0.2 0.5]});
%%
for i=2
    runid = sprintf('SN_Gbstar_%d_run', i);
    prob = coco_prob();
    
    prob = ode_ep2SN(prob, '', 'bdeps0244', SNlab(i));

    % add user points
    prob = coco_add_event(prob, 'UZ', 'b_star', 0.2:0.001:0.5);%linspace(0.1,0.5,20)
    coco(prob,runid , [], 1, {'G' 'b_star'}, {[0.3 4] [0.2 0.5]});
end
