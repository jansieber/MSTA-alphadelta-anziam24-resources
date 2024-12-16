dimpar.A= 3.2500;
dimpar.B= 22;
dimpar.C= 135;
dimpar.C1= 135;
dimpar.C2= 108;
dimpar.C3= 33.7500;
dimpar.C4= 33.7500;
dimpar.P= 0;
dimpar.a= 100;
dimpar.b= 50;
dimpar.e0= 2.5000;
dimpar.r= 0.5600;
dimpar.v0= 6;
% Example data
ydim = linspace(-10,20,300);
dimpar=load([pwd(),filesep,'dimeparam.mat']);% dowonload param scribt file 
sigm=@(x)2*dimpar.e0./(1+exp(dimpar.r.*(dimpar.v0-x)));
%%
bstar=dimpar.b/dimpar.a;
beta=2*dimpar.a/(dimpar.B*dimpar.r*dimpar.C*2*dimpar.e0);% 0.024;
epsilon=beta;
% epsilon=0.001
yndim= dimpar.r*dimpar.C*epsilon*ydim;%shall use beta?)
alpha1=1;
alpha2=0.8;
alpha3=0.25; alpha4=0.25;
c=dimpar.r*dimpar.v0;
% switch point
y01=c*beta;
y02=c*beta/alpha3;
y03=c*beta/alpha1;
a1=1;
a2=alpha3;
a3=alpha1;
% Calculate the sigmoid function
newsigm=@(x,x0,epsilon)1./(1+exp(((x0)-x)./(epsilon)));
%%
yrg=linspace(-0.1,0.4,1000);
ynrg=dimpar.r*dimpar.C*epsilon*yrg;
lw={'linewidth',2};
ltx={'Interpreter','latex'};
fna={'FontName','Courier New'};
fns={'FontSize',15};
fw={'FontWeight','bold'};
clr=lines();
figure(1);clf;
tl=tiledlayout(1,1);
ax1=axes(tl);
plot(ax1,dimpar.C1*yrg,sigm(dimpar.C1*yrg),'color',clr(1,:),lw{:});
xlabel(ax1, 'activation input $C_1 Y_1$ (mV)',ltx{:}); % X-label for ax2
ylabel(ax1, '$\mathrm{Sigm}(C_1\,Y)$',ltx{:}, 'Color', 'k'); % Y-label for ax2
xlim(ax1,dimpar.C1*yrg([1,end]))
%xl1=xline(ax1,dimpar.v0,'-',lw{:},'DisplayName','$v_0$','Color',clr(1,:));
%yl1=yline(ax1,dimpar.e0,'-',lw{:},'Color',(clr(1,:)+[1,1,1])/2,'DisplayName','$e_0$');
set(ax1,'LineWidth',1,fna{:},fns{:},fw{:},lw{:},...
    'XColor',clr(1,:),'YColor',clr(1,:),'box','off');
ax2=axes(tl);
plot(ax2,ynrg,newsigm(ynrg,y01,epsilon),'k--',lw{:});
xl2=xline(ax2,y01,'k:',lw{:},'DisplayName',sprintf('$v_0$ (dim.) $\\sim y_{0,1}$ (non-dim.)'));
yl2=yline(ax2,0.5,'-.',lw{:},'Color',0.5*[1,1,1],'DisplayName',sprintf('$e_0$ (dim.) $\\sim 1/2$ (non-dim.)'));
ax2.XAxisLocation='top';
ax2.YAxisLocation='right';
ax2.Color='none';
set(ax2,'LineWidth',1,fna{:},fns{:},fw{:},lw{:},...
    'XColor','k','YColor','k','box','off');
xlim(ax2,ynrg([1,end]))
ylim(ax2,[0,1])
xlabel(ax2, '$y_1$ (non-dim.)',ltx{:}); % X-label for ax2
ylabel(ax2, '$\mathrm{S}_{\epsilon/a_1}(y_1-y_{0,1})$',ltx{:}, 'Color', 'k'); % Y-label for ax2
legend(ax2,[xl2,yl2],'Location','SouthEast','FontSize',15,ltx{:})
%%
% exportgraphics(figure(1),'sigmoid.pdf')
