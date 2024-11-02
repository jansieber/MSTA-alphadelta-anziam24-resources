% xeq anaylsis for small epsilon=0.001
dimpar=load('dimeparam.mat');
beta=2*dimpar.a/(dimpar.B*dimpar.r*dimpar.C*2*dimpar.e0);% 0.024;
epsilon=0.001;
alpha1=1;
alpha2=0.8;
alpha3=0.25; 
alpha4=0.25;
c=dimpar.r*dimpar.v0;
% switch point
 y01=c*beta;
y02=c*beta/alpha3;
y03=c*beta;
a1=1;
a2=alpha3;
a3=alpha1;
b_star=dimpar.b/dimpar.a;
jnew= 2;

% state variable y1 
y1rag=linspace(0,3,5000000);
% dimensionless sigm function with epsilon
newsigm=@(y,y0,epsilon)1./(1+exp(((y0)-y)./(epsilon)));
%%
% G bifurcation paramter used in figure 3 

%  Hopf bifurcation  and SN is from data smalleps=001
 parnames={'jnew','G','P','alpha1','alpha2','alpha3','alpha4','b_star','c','beta',...
      'epsilon','x01','x02','x03','a1','a2','a3'};
 paramHB=bd_getpoint('smalleps=001','HB',parnames);
  paramSN=bd_getpoint('smalleps=001','SN',parnames);

%    G=paramHB(2,1);1.481791877233071%Hopf
%    G=paramHB(2,2);
%    G=paramSN(2);

G_values =[1.4  1.481791877233071    3.5  paramHB(2,2)    9   paramSN(2)      22];

%G=jnew*alpha2/y01;% SN1 
%G=jnew/y03;
%G=jnew/y02; %Hopf2
%G=b_star*alpha2*jnew/(b_star*y01+(alpha4*jnew));% Hopf1 

n_tiles = length(G_values);
tiledlayout(1, n_tiles+1); % Adjust rows/columns as needed

% loop to go over each value 

for i = 1:n_tiles
    G = G_values(i);

% xeq_y2 expression

y2=@(x)jnew*alpha4/(b_star)*newsigm(x,y02,epsilon/a2);

% xeq_y3 expression

y3=@(x1)(jnew*alpha2/(G))*newsigm(x1,y03,epsilon/a3);%S(..,y03,.epsilon)

y=@(x) y3(x)- y2(x);
% r.h.s of Y1
f=@(x)newsigm(y(x),y01,epsilon/a1);

% L.h.s of Y1
l = @(x) (G/jnew)*x;
nexttile;

fna={'FontName','Courier New'};
fns={'FontSize',14};
fw={'FontWeight',...
    'bold'};
ltx={'Interpreter','latex'};
set(gca,'LineWidth',1,fna{:},fns{:},fw{:}...
    );

plot(y1rag,f(y1rag),'LineWidth', 2);
hold on;grid on;
plot(y1rag,l(y1rag), 'LineWidth', 2);

plot([y03, y03], [0, 1.2], 'k--', 'LineWidth', 1.5); 
plot([y02, y02], [0, 1.2], 'k--', 'LineWidth', 1.5); 
plot(y03, 0, 'g^', 'MarkerSize', 12, 'MarkerFaceColor', 'g'); 
plot(y02, 0,'^', 'MarkerSize', 12,'color',[0.49,0.18,0.56],'MarkerFaceColor', [0.49,0.18,0.56]);

title(['${G=' num2str(G) '}$'], 'Interpreter', 'latex');

xlim([0, 1.47]);            
ylim([0 1.2]);


% intersection point 
% Define RHS - LHS
h = @(x) f(x) - l(x);
% Find the intersection point (initial guess at x = 0)
x_intersect = fzero(h, 0);
y_intersect = f(x_intersect);
hold on;
plot(x_intersect, y_intersect, 'ko', 'MarkerSize',10, 'MarkerFaceColor','k');

x_intersect1 = fzero(h, y03);
y_intersect1 = f(x_intersect1);

plot(x_intersect1, y_intersect1, 'ko', 'MarkerSize',10, 'MarkerEdgeColor','k'); 

x_intersect2 = fzero(h, y02);
y_intersect2 = f(x_intersect2);
plot(x_intersect2, y_intersect2, 'ko', 'MarkerSize',10, 'MarkerEdgeColor','k'); 
x_intersect3 = fzero(h, jnew/G);

y_intersect3 = f(x_intersect3);
plot(x_intersect3, y_intersect3, 'ko', 'MarkerSize',10, 'MarkerFaceColor','k');
end

%%
legend('$\mathrm{S_{\epsilon/a_1}((y_3-y_2)-y_{01})}$','${(G/2)y_1}$','','','${y_{03}}$','${y_{02}}$','Stable','Unstable','interpreter','latex');
%%
% plot one para bifurcation in (y1,G) plane 
nexttile

epplot=ep_plot_theme('ep');
epplot = struct('special', {{'SN', 'HB'}});
epplot.SN={'kd' , 'MarkerFaceColor' , 'g'  ,'MarkerSize' , 14};
epplot.HB={'kd' , 'MarkerFaceColor' , 'r'  ,'MarkerSize' , 14};
epplot.lspec{1}={'b-' ,   'LineWidth',   2};
epplot.lspec{2}={'b--'  ,  'LineWidth' ,   2};
coco_plot_bd(epplot,'smalleps=001' , 'G','x',@(x)x(1,:));
ylabel('${y_1}$(pc)', 'Interpreter', 'latex');
xlabel('${G}$', 'Interpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

% Set the x-ticks 
xticks(ax, [0, 1.482080013460666, 3, 5, 6.325289511885655,18, 19.405212300224832]);

% Set the corresponding x-tick labels
xticklabels(ax, {'$0$', '$G_1$', '$3$', '$5$', '$G_2$',     '18',  '$G_3$'});
yticks(ax, [0, x01, x02, 1, 1.4]);
yticklabels(ax, {'$0$', '$y01$', '$y02$', '$1$', '$1.4$'});
set(gca, 'TickLabelInterpreter', 'latex');
%%
 saveas(figure(1),'fig3.svg', 'svg');
