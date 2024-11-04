figure;

%%
% plot colormap 

nexttile;                        % Move to the first tile
ax=gca;
s=load('Ghopf.mat');          Ghopf      =s.Ghopf;
s=load('Gsn.mat');            Gsn        =s.Gsn;
s=load('G_vals6050.mat');     G_vals     =s.G_vals;
s=load('bs_vals6050.mat');    bs_vals    =s.bs_vals;
s=load('period_vals6050.mat');period_vals=(s.period_vals);
a=100; 
fre= a*(1./period_vals);
clr=colormap('summer');
clr=clr(end:-1:1,:);
cvals=linspace(0,1,size(clr,1));
pp=interp1(cvals,clr,'linear','pp');
clr_res=ppval(pp,1-cvals.^(1/2));
[cmap,cf]=contourf(ax,bs_vals,G_vals,fre,200,'Linestyle','none');
ax.Colormap=clr_res;
ax.CLim(1)=0;% lower bound 
cb=colorbar(ax);
cb.Ticks=round(ax.CLim(1)):2:round(ax.CLim(2),-1);
cb.Label.String='Frequency $(Hz)$';
cb.Label.Interpreter='latex';
cb.Label.FontSize=18;
hold on;
phopf=plot(ax,bs_vals(1,:),Ghopf,'r','DisplayName','Hopf bifurcation','LineWidth',2);
psn=plot(ax,bs_vals(1,:),Gsn,'b','DisplayName','Saddle node (SNIC)','LineWidth',2);
ltx={'Interpreter','latex'};
legend(ax,[phopf,psn],'Location','best');
set(ax,'Fontsize',25,'FontName','courier','FontWeight','bold','LineWidth',1)

%%
%% plot grazing bifurcation curve 

hold on; 
varnames={'y01','y02','y03','bstar','jnew','G','alpha2','alpha4',...
    'ts1off','ts1on','ts2off','ts2on','T',...
    't1min','y1_min','y1dot_min'};
cl=[varnames;num2cell(1:length(varnames))];
iv=struct(cl{:});                % automatically number parameters
filename='bstar_G';

bG=coco_bd_read(filename);   % read from data/bd_p/bd.mat
ubg=coco_bd_col(bG,varnames);

graze=plot(ax,ubg(iv.bstar,:),ubg(iv.G,:),'DisplayName','Grazing bifurcation','color','k','LineWidth',2);
%%
% extraxt grazing from data bstar_curve
filename='bstar_curve';
bda=coco_bd_read(filename);
u01=coco_bd_col(bda,varnames);
cl=[varnames;num2cell(1:length(varnames))];
iv=struct(cl{:}); 
upar=bd_getpoint('bstar_curve','graze',varnames);
plot(upar(iv.bstar),upar(iv.G),'kp', 'MarkerSize',14, 'MarkerFaceColor','g');
%%
% sn of periodic orbits 
hold on; 
for i=1:2
thm = ep_plot_theme('ep.SN');
thm.lspec={'m',  'LineWidth',  2};
coco_plot_bd(thm, sprintf('Posn_%d_run',i),'b_star', 'G');
end
xlabel(ax,'$b^*$',ltx{:});
ylabel(ax,'$G$',ltx{:});
xlim([0.2 0.5]);

%% sn of xeq for different epsilon 
varnames={'bstar','jnew','G','alpha2','alpha4',...
    'a1','a2','a3','y01','y02','y03','y1','epsilon'...
    'eqb','eqbd'};
format compact
cl=[varnames;num2cell(1:length(varnames))];
iv=struct(cl{:});
filename='sn_Gbstareps015';
filename1='sn_Gbstareps01';
%bG=coco_bd_read(filename);
bG1=coco_bd_read(filename);
bG2=coco_bd_read(filename1);
ubg1=coco_bd_col(bG1,varnames);
ubg2=coco_bd_col(bG2,varnames);
sn1=plot(ax,ubg1(iv.bstar,:),ubg1(iv.G,:),'blu--','DisplayName','SNIC,$\epsilon=0.015$','LineWidth',2);
sn2=plot(ax,ubg2(iv.bstar,:),ubg2(iv.G,:),':blu','DisplayName','SNIC,$\epsilon=0.01$','LineWidth',2);
%% time profile 

nexttile;    
ax=gca;
% time profile of alpha 
parnames={'jnew','G','P','alpha1','alpha2','alpha3','alpha4','b_star','c','beta',...
    'epsilon','x01','x02','x03','a1','a2','a3'};
format compact
cl=[parnames;num2cell(1:length(parnames))];
ip=struct(cl{:});               
F=sco_gen(@symnondritjan_eps);       
funcs = {F(''),F('x'),F('p')};   % r.h.s, derivative rt x and derivative wrt p
% dimpar=load('dimeparam.mat');
beta=2*dimpar.a/(dimpar.B*dimpar.r*dimpar.C*2*dimpar.e0);% 0.024;
alpha1=1; alpha3=0.25;
jnew=2;
c=dimpar.r*dimpar.v0;
% here call y01=x01
% y02=x02
% y03=y03
x01=c*beta;
x02=c*beta/alpha3;
x03=c*beta/alpha1;
a1=1;
a2=alpha3;
a3=alpha1;
par([ip.x01;ip.x02;ip.x03;ip.a1;ip.a2;ip.a3;ip.beta;ip.epsilon])=...
    [x01;     x02;    x03; a1;   a2;    a3;   beta;  0.024];
par([ip.jnew; ip.P;ip.alpha1;ip.alpha2;ip.alpha3;ip.alpha4;ip.c ])=...
    [jnew;       0;        alpha1; 0.8;  alpha3;    0.25; c];
dim=7;
%x0=zeros(dim,1);
x0=[0.9243;  1.1141;  1.1428; -0.3743;  0.0369;  0.0001;  0.0287];
G=1.4;
bstar=0.4;
par([ip.G;ip.b_star])=[G;bstar];
u0=par';
f_ode=@(t,x)funcs{1}(x,u0);
tspan=[0,100];
sol=ode45(f_ode,tspan,x0,odeset('RelTol',1e-6));
% change timestep
t=linspace(0,tspan(end),round((tspan(end)-tspan(1))/0.01));

lw_style = {'LineWidth', 2, 'LineStyle', '--'};
ltx={'Interpreter','latex'};
 hold on;
PC=plot(ax,t,deval(sol,t,1),'DisplayName','$y_1$(pc)','Color',[0.30,0.75,0.93],'LineWidth',2);
IN=plot(ax,t,(deval(sol,t,2)),'DisplayName','$y_2(inh)$','Color',[0.85,0.33,0.10],'LineWidth',2);
Ex=plot(ax,t,(deval(sol,t,3)),'DisplayName','$y_3$(exc)','Color',[0.47,0.67,0.19],'LineWidth',2);
yline(x03, lw_style{:});
legend(ax,[PC,IN,Ex],'Location','best');
set(ax,'Fontsize',25,'FontName','courier','FontWeight','bold','LineWidth',1)
xlabel('Time',ltx{:});
text(23,1.3, ' $\alpha$-like wave ', 'FontSize',25,ltx{:});
text(20,x03, ' $y_{03}$', 'FontSize',25,ltx{:})
xlim([0,30]);
title('a');
%% delta activity
nexttile
G=1.6;
bstar=0.4;
%x0=zeros(dim,1);
x0=[0.9243;  1.1141;  1.1428; -0.3743;  0.0369;  0.0001;  0.0287];
par([ip.G;ip.b_star])=[G;bstar];
u0=par';
f_ode=@(t,x)funcs{1}(x,u0);
tspan=[0,100];
sol=ode45(f_ode,tspan,x0,odeset('RelTol',1e-6));

PC=plot(ax,t,deval(sol,t,1),'DisplayName','$y_1$(pc)','Color',[0.30,0.75,0.93],'LineWidth',2);
IN=plot(ax,t,(deval(sol,t,2)),'DisplayName','$y_2(inh)$','Color',[0.85,0.33,0.10],'LineWidth',2);
Ex=plot(ax,t,(deval(sol,t,3)),'DisplayName','$y_3$(exc)','Color',[0.47,0.67,0.19],'LineWidth',2);
yline(x03, lw_style{:});

legend(ax,[PC,IN,Ex],'Location','best');
set(ax,'Fontsize',25,'FontName','courier','FontWeight','bold','LineWidth',1)
ylabel(ax,'Neural activity',ltx{:});
xlabel(ax,'Time',ltx{:});
text(23,1.3, ' $\delta$-like wave', 'FontSize',25,ltx{:});
% text(0.4,1.6, ' $\delta$', 'FontSize',18, 'FontWeight','bold',ltx{:})
xlim([0,30]);
%%
