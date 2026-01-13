% time profile of alpha canard Hopf bifurcation 

% take point near Hopf bifurcation G1
fileread=coco_bd_read('po_Canardeps001UZ');
labs=coco_bd_labs(fileread,'uz');
profile=po_read_solution('po_Canardeps001UZ',labs(81)); % cANRAD 
xnew=[profile.xbp(1:end-1,:);profile.xbp];
tnew=[profile.tbp(1:end-1);profile.tbp+profile.T];
lw={'linewidth',2};

%%

%%
dimpar=load([pwd(),filesep,'dimeparam.mat']);% dowonload param scribt file ;
alpha1=1;alpha2=0.8; alpha3=0.25;alpha4=0.25;
beta=2*dimpar.a/(dimpar.B*dimpar.r*dimpar.C*2*dimpar.e0);% 0.024;
epsilon=0.001;
jnew=2;
c=dimpar.r*dimpar.v0;
y01=c*beta;
y02=c*beta/alpha3;                                                                                                                                                          
y03=c*beta/alpha1;
a1=1;
a2=alpha3;
a3=alpha1;
% sigm function with epsilon
newsigm=@(y,y0,epsilon)1./(1+exp(((y0)-y)./(epsilon)));
%%

figure;
t=tiledlayout(1,2);
ax1=nexttile(1);
ltx={'Interpreter','latex'};
fna={'FontName','Courier New'};
fns={'FontSize',18};
fw={'FontWeight','bold'};
set(gca,'LineWidth',1,fna{:},fns{:},fw{:});
box(ax1, 'on');

x1=plot(ax1,tnew,xnew(:,1),lw{:},'Color',[0.30,0.75,0.93],'DisplayName','$y_{1}$(pc)');
hold on; 
x2=plot(tnew,newsigm(xnew(:,3)-xnew(:,2),y01,(epsilon/a1)),lw{:},'DisplayName','$\mathrm{S}_\epsilon/a_1$');
x2. LineStyle={'--'};
% y03
x3=yline(ax1,y03, 'LineWidth', 2,'Color',[0.50,0.50,0.50], 'DisplayName','$y_{0,3}$'); %  horizontal line at y = y03

x4=yline(ax1,y02, 'LineWidth', 2,'Color',[0.93,0.69,0.13], 'DisplayName','$y_{0,2}$'); % horizontal line at y = y02

x5=plot(ax1,tnew,xnew(:,3),'color',[0.47,0.67,0.19],lw{:},'DisplayName','$y_{3}$(exc)');

legend(ax1,[x1,x2,x3,x4,x5],'Location','EastOutside','FontSize',20,ltx{:})
xlabel('$t$','Interpreter','latex');
xlim([0,profile.T]);

%%
ax2=nexttile(2);
% y2
z1=plot(ax2,tnew,xnew(:,2),lw{:},'Color',[0.85,0.33,0.10],'DisplayName','$y_{2}$(inh)');
set(gca,'Fontsize',14,'FontName','courier','FontWeight','bold','LineWidth',1)

hold on; 
z2=plot(ax2,tnew,newsigm(xnew(:,1),y02,(epsilon/a2)),'--',lw{:},'Color',[0.64,0.08,0.18],'DisplayName','$\mathrm{S}_\epsilon/a_2$');

y1thr=xnew(1,3)-y01;

z3=yline(ax2,y1thr,lw{:},'Color',[1.00,0.07,0.65],'DisplayName','$y_3-y_{0,1}$'); % Add a horizontal line at y = y3-y01

legend(ax2,[z1,z2,z3],'Location','EastOutside','FontSize',20,ltx{:})
legend('$y_{2}$(inh)','$\mathrm{S}_\epsilon/a_2$','$y_3-y_{0,1}$',ltx{:},'Location','EastOutside','FontSize',20,ltx{:})
xlim([0,profile.T]);

% t.TileSpacing = 'compact'; 
% t.Padding = 'compact';  
%%
saveas(figure(5),'alphalikecanardupdata2.svg','svg');
