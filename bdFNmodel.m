clear
parnames={'a','epsilon'};
format compact
cl=[parnames;num2cell(1:length(parnames))];
ip=struct(cl{:});                % automaticaly number parameters
F=sco_gen(@sym_FN);       % F and F2 are the same
funcs = {F(''),F('x'),F('p')};   % r.h.s, derivative rt x and derivative wrt p
a=1.2; epsilon=0.2;
par([ip.a;ip.epsilon])=...
    [a;  epsilon];
dim=2;
x0=zeros(dim,1);
prob = coco_prob();
prob = ode_isol2ep(prob, '', funcs{:}, x0, parnames, par);
prob = coco_set(prob, 'cont', 'PtMX', 5000,'h_max',0.01,'h_min',1e-6,'NPR',100);

fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
    'ep_run');
% start from a=1.2 then up to 2
% again start from a=1.2 then down  to -1
 coco(prob, 'bdFN', [], 1, 'a', [-1 2]);

%%

%% branch off at Hopf
bd_ep_P_0=coco_bd_read('bdFN');
HBlab = coco_bd_labs(bd_ep_P_0, 'HB');
prob = coco_prob();
prob = coco_set(prob, 'coll', 'MXCL',false); % remove discretization error limit
prob = ode_HB2po(prob, '', 'bdFN', HBlab(2));
prob = coco_add_event(prob, 'UZ', 'a',-1:0.1:0.5);%-1.0000  -0.9000...-0.5000
prob = coco_set(prob, 'cont', 'PtMX',[700,0], 'NAdapt', 5,...
   'h_max',1e3,'NPR',10,'MaxRes',100,'norm',inf);
coco(prob, 'po_HB1', [], 1, {'a','po.period'}, {[-1 2],[0 510]});
% this data used to produce time profile in Fig 2

%% copy here 
for i=2
    runid = sprintf('SN_Gbstar_%d_run', i);
    prob = coco_prob();
    
    prob = ode_ep2SN(prob, '', 'bdeps0244', SNlab(i));

    % add user points
    prob = coco_add_event(prob, 'UZ', 'b_star', 0.2:0.001:0.5);%linspace(0.1,0.5,20)
    coco(prob,runid , [], 1, {'G' 'b_star'}, {[0.3 4] [0.2 0.5]});
end



%%

copy here 
%% Start continuation along family of Hopf bifurcations

% parameters 'a2' and 'G'
% read data of ep solution for varying G 
bd_a2=coco_bd_read('epa2=alph3');
% get Hopf points 
HBlab = coco_bd_labs(bd_a2, 'HB');

% loop over each Hopf on the upper branch of Hopf 
for i=HBlab(2)
  runid = sprintf('HBG_a2_%d_run', i);
  prob = coco_prob();
 
  prob = ode_ep2HB(prob, '', 'epa2=alph3', i); 
  
  % add user points
  %prob = coco_add_event(prob, 'UZ', 'a2', linspace(0.01,0.65,20));

   prob = coco_add_event(prob, 'UZ', 'a2', linspace(0.25,0.5,6));
   prob = coco_add_event(prob, 'Takens', 'boundary','ep.test.BTP','==', 0);
   fprintf(...
  '\n Run=''%s'': Continue family of periodic orbits from point %d in run ''%s''.\n', ...
  runid, i, 'HB_run');

  coco(prob, runid, [], 1,  {'a2' 'G'}, {[0.01 1] [0.5 40]});
end

%%

clr=lines();
for i=1:60
labs=coco_bd_labs('po_HB1','UZ');
profile=po_read_solution('po_HB1',labs(i));
t=profile.tbp;
profile.p
figure
plot(profile.xbp(:,1),...
    profile.xbp(:,2),'Color',clr(2,:));
end 
%%
plot(t,profile.xbp(:,1),'Color',clr(2,:),...
    t,profile.xbp(:,2),'Color',clr(2,:)); 
figure
plot(t,profile.xbp(:,1),'Color',clr(2,:))


%% phase plane 
plot(profile.xbp(:,1),...
    profile.xbp(:,2),'Color',clr(2,:));

%%
%%
hold on; 
thm=po_plot_theme('po');
thm = struct('special', {{'SN'}});
thm.SN={'kd' , 'MarkerFaceColor' , 'g'  ,'MarkerSize' , 14};
thm.lspec{1}={'k-' ,   'LineWidth',   2};
thm.lspec{2}={'k--'  ,  'LineWidth' ,   2};
coco_plot_bd(thm, 'po_HB1', 'a','MAX(x)',@(x)x(1,:))
coco_plot_bd(thm, 'po_HB1', 'a','MIN(x)',@(x)x(1,:))
% ylabel('$\mathrm{y_1}$(pc)', 'Interpreter', 'latex');
lw = {'LineWidth', 2};
ltx={'Interpreter','latex'};
fna={'FontName','Courier New'};
fns={'FontSize',25};
% Set the remaining axes properties
set(gca,'FontSize',14,'LineWidth',...
    1);
%%
%%  plot nullclines 

% Define the range for x
x = -2:0.1:2;  % from -10 to 10 with step size 0.1

% Compute y = x^3/3-x
y = (x.^3/3)-x;  % use element-wise power operator .^

% Plot the function
figure
plot(x, y);
xlabel('x');
ylabel('y = x^3');
title('Plot of y = x^3');
grid on;