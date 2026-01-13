% run by coco_March 2020
%%
% run startup function of coco_March 2020
% run file gen_sys_nodimwitheps

%% Load path
% clear
%cocofolder=fileparts(which('coco'));
%symcocopath=[cocofolder,'/../../contributed/symcoco/toolbox'];
%addpath(symcocopath); % path of symcoco routines
parnames={'jnew','G','P','alpha1','alpha2','alpha3','alpha4','b_star','c','beta',...
    'epsilon','x01','x02','x03','a1','a2','a3'};
format compact
cl=[parnames;num2cell(1:length(parnames))];
ip=struct(cl{:});                % automaticaly number parameters
F=sco_gen(@symnondritjan_eps);       % F and F2 are the same
funcs = {F(''),F('x'),F('p')};   % r.h.s, derivative rt x and derivative wrt p
%% initial guess for parameters
dimpar=load('dimeparam.mat');
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
prob = ode_isol2ep(prob, '', funcs{:}, x0, parnames, double(par));
prob = coco_set(prob, 'cont', 'PtMX',  [50  60000],'h0',1e-4,'h_min',1e-6,'h_max',0.05,'NPR',500);
fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
  'ep_run');
coco(prob, 'smalleps=001', [], 1, 'G', [0.3 20]);


%% max and min
bd_ep_P_0=coco_bd_read('smalleps=001');
HBlab = coco_bd_labs(bd_ep_P_0, 'HB');
prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST',100,'NTSTMX',200); % remove discretization error limit
%prob = coco_set(prob, 'po', 'bifus',false); % disable bif detection
prob = ode_HB2po(prob, '', 'smalleps=001', HBlab(1));
prob = coco_set(prob, 'cont', 'PtMX',[0,500], 'NAdapt', 1,...
   'h_max',1e5,'NPR',10,'MaxRes',100,'norm',inf);
prob = coco_add_event(prob, 'uz', 'G', 1.4821:0.0001:1.5);
coco(prob, 'po_Canardeps001UZ', [], 1, {'G','po.period'}, {[0.5 10],[0,50]});
%coco(prob, 'po_smalleps001UZ', [], 1, 'G', [0.5 4]);
%% deco for alpha oscillations and grazing
%% max and min
bd_ep_P_0=coco_bd_read('smalleps=001');
HBlab = coco_bd_labs(bd_ep_P_0, 'HB');
prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST',100,'NTSTMX',200); % remove discretization error limit
%prob = coco_set(prob, 'po', 'bifus',false); % disable bif detection
prob = ode_HB2po(prob, '', 'smalleps=001', HBlab(1));
prob = coco_set(prob, 'cont', 'PtMX',[0,500], 'NAdapt', 1,...
   'h_max',1e5,'NPR',10,'MaxRes',100,'norm',inf);
prob = coco_add_event(prob, 'UZ', 'G', [1.7,2.05,2.1]);
coco(prob, 'po_grazing_illu', [], 1, {'G','po.period'}, {[0.5 2.101],[0,50]});
%%
bd=coco_bd_table('po_grazing_illu','numlab',true);
mx=cat(2,bd.('MAX(x)'){:});
mn=cat(2,bd.('MIN(x)'){:});
clr=lines();
axprops={'LineWidth',1,'box','on','FontSize',16};
lw={'LineWidth',2};
g_rg=[1.4,2.12];
fig=figure(1);clf;tiledlayout(2,2,'TileSpacing','tight');
ax1=nexttile(1);ax2=nexttile(3);
%figure(2);clf;ax=gca;
hold(ax1,'on');hold(ax2,'on');
thm1=po_plot_theme('po');
thm3=thm;
thm1.lspec{1}([1,4:7])={'o-','MarkerSize',4,'Color',clr(1,:)};
thm3.lspec{1}([1,4,5])={'o-','Color',clr(2,:)};
coco_plot_bd(ax1,thm1,'po_grazing_illu','G','MAX(x)',1);
coco_plot_bd(ax1,thm1,'po_grazing_illu','G','MIN(x)',1);
coco_plot_bd(ax2,thm3,'po_grazing_illu','G','MAX(x)',3);
coco_plot_bd(ax2,thm3,'po_grazing_illu','G','MIN(x)',3);
thme1=struct('special',{{'SN','HB'}});
thme3=thme1;
thme3.lspec{1}={'-','Color',clr(2,:)};
thme3.lspec{2}={'--','Color',clr(2,:)};
coco_plot_bd(ax1,thme1,'smalleps=001','G','MAX(x)',1);
coco_plot_bd(ax2,thme3,'smalleps=001','G','MAX(x)',3);
yline(ax1,bd.x03(1),'k:','LineWidth',2);
i_gr=find(mn(1,:)<bd.x03(1),1,'first')-2;
plot(ax1,bd.G(i_gr),mn(1,i_gr),'ko','MarkerFaceColor','k');
uzprop={'-.','color',clr(4,:),'LineWidth',1};
xline(ax1,bd.G(i_gr),uzprop{:});
xline(ax2,bd.G(i_gr),uzprop{:});
xlim(ax1,g_rg);
ax1.XTickLabel={};
xlim(ax2,g_rg)
i_uz=find(strcmp(bd.TYPE,'UZ'));
xline(ax1,bd.G(i_uz([1,3])),uzprop{:})
xline(ax2,bd.G(i_uz([1,3])),uzprop{:})
text(ax1,g_rg(1)+0.01,bd.x03(1),'$y_{0,3}$','Interpreter','latex',...
    'FontSize',18,'VerticalAlignment','bottom');
ylabel(ax1,'$[\min y_1,\max y_1]$','Interpreter','latex')
ylabel(ax2,'$[\min y_3,\max y_3]$','Interpreter','latex')
xlabel(ax2,'$G$','Interpreter','latex')
set(ax1,axprops{:});
set(ax2,axprops{:});
hopf=findobj(ax2,'DisplayName','Hopf Bifurcation');
ax3=nexttile(2);hold(ax3,'on');
ax4=nexttile(4);hold(ax4,'on');
rg=bd.LAB(find(strcmp(bd.TYPE,'UZ')));
clear po
for i=1:length(rg)
    po(i)=po_read_solution('','po_grazing_illu',rg(i));
end
for i=1:length(rg)
    fac=1/1.1^i;
    y1=plot(ax3,po(i).tbp/po(i).T,po(i).xbp(:,1),'-','Color',clr(1,:)*fac,lw{:});
    y3=plot(ax4,po(i).tbp/po(i).T,po(i).xbp(:,3),'-','Color',clr(2,:)*fac,lw{:});
    u3=plot(ax3,po(i).tbp/po(i).T,newsigm(po(i).xbp(:,1),po(i).p(ip.x03),po(i).p(ip.epsilon)),...
        '-','Color',clr(3,:)*fac,lw{:});
end
yline(ax3,bd.x03(1),'k:','LineWidth',2);
%ylabel(ax3,'$y_1(t)$','Interpreter','latex')
%ylabel(ax4,'$y_3(t)$','Interpreter','latex');
xlabel(ax4,'time $t/T$','Interpreter','latex')
set(ax3,axprops{:});
set(ax4,axprops{:});
ax3.XTickLabel={};
ylim(ax4,[0.6,1]);
ylim(ax2,[0.6,1.2]);
ylim(ax3,ax1.YLim);
ax3.YTick=ax3.YTick(1:end-1);
ax3.YTickLabel=cellfun(@(c){['   ',c]},ax3.YTickLabel);
text(ax3,0.5,bd.x03(1),'$y_{0,3}$','Interpreter','latex',...
    'FontSize',18,'VerticalAlignment','bottom');
lg=legend(ax3,[hopf,y1,y3,u3],{'Hopf','$y_1(t)$','$y_3(t)$','$u_3(t)$'},'Interpreter','latex');
lg.Layout.Tile='east';
rel=@(bd,r,i)bd(i)+(bd(2)-bd(1))*r;
labprops={'Interpreter','latex','FontSize',18};
text(ax1,rel(ax1.XLim,-0.28,1),rel(ax1.YLim,0,2),'a)',labprops{:})
text(ax2,rel(ax2.XLim,-0.28,1),rel(ax2.YLim,0,2),'c)',labprops{:})
text(ax3,rel(ax3.XLim,-0.12,1),rel(ax3.YLim,0,2),'b)',labprops{:})
text(ax4,rel(ax4.XLim,-0.12,1),rel(ax4.YLim,0,2),'d)',labprops{:})
get(fig,'Position');
fig.Position(3:4)=[775,340];
%%
exportgraphics(fig,'../PhD-Thesis-draft/figures/po_grazing_illu.pdf','ContentType','vector');