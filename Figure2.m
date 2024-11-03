%% Plot bifurcation diagram as function of A
%%
% Y3 versus  A
dimfile='ep_inputp0'

figure(1);clf; hold on;
grid on; box on;
lw = {'LineWidth', 2};
fna={'FontName','Courier New'};
fns={'FontSize',14};
fw={'FontWeight',...
    'bold'};
ltx={'Interpreter','latex'};
set(gca,'LineWidth',1,fna{:},fns{:},fw{:}...
    );
tiledlayout(3,3);
nexttile(1);
epplot=ep_plot_theme('ep');
epplot = struct('special', {{'SN', 'HB'}});
epplot.SN={'kd' , 'MarkerFaceColor' , 'g'  ,'MarkerSize' , 14};
epplot.HB={'kd' , 'MarkerFaceColor' , 'r'  ,'MarkerSize' , 14};
epplot.lspec{1}={'b-' ,   'LineWidth',   2};
epplot.lspec{2}={'b--'  ,  'LineWidth' ,   2};
coco_plot_bd(epplot, dimfile, 'A','x',@(x)x(3,:));

hold on;
thm=po_plot_theme('po');
thm = struct('special', {{'SN'}});
thm.SN={'kd' , 'MarkerFaceColor' , 'g'  ,'MarkerSize' , 14};
thm.lspec{1}={'k-' ,   'LineWidth',   2};
thm.lspec{2}={'k--'  ,  'LineWidth' ,   2};
% periodic  orbit
dimpo='po_inputp0';
coco_plot_bd(thm, dimpo, 'A','MAX(x)',@(x)x(3,:));
coco_plot_bd(thm, dimpo, 'A','MIN(x)',@(x)x(3,:));
ylabel('${Y_3}$(exc)', 'Interpreter', 'latex');
%%
% freq vs A
% read data
dimpo='po_inputp0';
bd_po=coco_bd_read(dimpo);
parnames={'A','B','a','b','C1','e0','r','v0','P','po.period'};
ubg=coco_bd_col(bd_po,parnames);
%figure; 
nexttile(2);
hold on;
plot(ubg(1,:),1./ubg(10,:),lw{:});
freq=1./ubg(10,:);
%  convert the bd's to tables and use semilogy
bd01=cell2table(bd_po(2:end,:),'VariableNames',bd_po(1,:));
% update colunm of period (replace column values with freq. value)
bd01.('po.period')=freq';
geti=@(x,i)x(i);
indices=[1;find(strcmp(bd01.('TYPE'),'SN'));length(bd01.('TYPE'))];
% nexttile(2);
for i=1:length(indices)-1
    rg=indices(i):indices(i+1);
    if i==1
        plot(geti(bd01.('A'),rg),geti(bd01.('po.period'),rg),'k','LineWidth',2);
    elseif i==2
        plot(geti(bd01.('A'),rg),geti(bd01.('po.period'),rg),'k--','LineWidth',2);
    elseif i==3
        plot(geti(bd01.('A'),rg),geti(bd01.('po.period'),rg),'k','LineWidth',2);
        hold on;grid on; box on;
    end
end
% This should  give the indices of SN.
%%
indices=find(strcmp(bd01.('TYPE'),'SN'))
plot(geti(bd01.('A'),indices(1)),geti(bd01.('po.period'),indices(1)),'kd','MarkerFaceColor','g','MarkerSize',14);
plot(geti(bd01.('A'),indices(2)),geti(bd01.('po.period'),indices(2)),'kd','MarkerFaceColor','g','MarkerSize',14);

% ylabel('$\mathrm{Frequency\ (Hz)}$', 'Interpreter', 'latex');
% xlabel('${A}$', 'Interpreter', 'latex');

%%
% freq vs A for p=120
dimpo='poinputp120';
bd_po=coco_bd_read(dimpo);
parnames={'A','B','a','b','C1','e0','r','v0','P','po.period'};
ubg=coco_bd_col(bd_po,parnames);
hold on; 
plot(ubg(1,:),1./ubg(10,:),lw{:});
freq=1./ubg(10,:);
%  convert the bd's to tables and use semilogy
bd01=cell2table(bd_po(2:end,:),'VariableNames',bd_po(1,:));
% update colunm of period (replace column values with freq. value)
bd01.('po.period')=freq';
geti=@(x,i)x(i);
indices=[1;find(strcmp(bd01.('TYPE'),'SN'));length(bd01.('TYPE'))];
 nexttile(3);
for i=1:length(indices)-1
    rg=indices(i):indices(i+1);
    if i==1
        plot(geti(bd01.('A'),rg),geti(bd01.('po.period'),rg),'k','LineWidth',2);
    elseif i==2
        plot(geti(bd01.('A'),rg),geti(bd01.('po.period'),rg),'k--','LineWidth',2);
    elseif i==3
        plot(geti(bd01.('A'),rg),geti(bd01.('po.period'),rg),'k','LineWidth',2);
        hold on;grid on; box on;
    end
end
ylabel('$\mathrm{Frequency\ (Hz)}$', 'Interpreter', 'latex');
xlabel('${A}$', 'Interpreter', 'latex');
%%

%time profile for alpha && delta

%% Example for getting the solution profile (here alpha and delta activity )
% alpha, A=11
% take label point from data
profile=po_read_solution('po_inputp0',5);
% two y axis
nexttile;
[AX,H1,H2]=plotyy(profile.tbp,profile.xbp(:,2:3),profile.tbp,profile.xbp(:,1));
set(...
    get(AX(1),'Ylabel'),...
    'String', 'Neural activity',...
    'Fontsize',14 ...
    )
set(AX(1),'LineWidth',1,fna{:},fns{:},fw{:}...
    );
set(AX(1),'Xlim',[0 profile.T]);
set(AX(2),'Xlim',[0 profile.T]);


set(AX(1), 'YTick', [0 70], 'YTickLabel', {'0', '70'})
set(AX(2), 'YTick', [0 0.5], 'YTickLabel', {'0', '0.5'})
set(AX(2), 'Ylim', [0 0.5]);
set(AX(1), 'Ylim', [0 70]);

set(H1,"LineWidth",2);
set(H1(1),'Color',[0.85,0.33,0.10])

set(H1(2),'Color',[0.47,0.67,0.19])

set(H2,'Color',[0.30,0.75,0.93])
set(H2,"LineWidth",2)

ltx={'Interpreter','latex'};
set(AX(2),'LineWidth',1,fna{:},fns{:},fw{:}...
    );
set(AX(2), 'YColor', [0.30,0.75,0.93]);
xlabel('$\mathrm{Time\ (sec)}$', 'Interpreter', 'latex');

%%
% delta, A=10
profile=po_read_solution('po_inputp0',17);

% two y axis
nexttile;
[AX,H1,H2]=plotyy(profile.tbp,profile.xbp(:,2:3),profile.tbp,profile.xbp(:,1));
set(...
    get(AX(1),'Ylabel'),...
    'String', 'Neural activity',...
    'Fontsize',14 ...
    )
set(AX(1),'LineWidth',1,fna{:},fns{:},fw{:}...
    );
set(AX(1),'Xlim',[0 profile.T]);
set(AX(2),'Xlim',[0 profile.T]);


set(AX(1), 'YTick', [0 70], 'YTickLabel', {'0', '70'})
set(AX(2), 'YTick', [0 0.5], 'YTickLabel', {'0', '0.5'})
set(AX(2), 'Ylim', [0 0.5]);
set(AX(1), 'Ylim', [0 70]);

set(H1,"LineWidth",2);
set(H1(1),'Color',[0.85,0.33,0.10])

set(H1(2),'Color',[0.47,0.67,0.19])

set(H2,'Color',[0.30,0.75,0.93])
set(H2,"LineWidth",2)

ltx={'Interpreter','latex'};
set(AX(2),'LineWidth',1,fna{:},fns{:},fw{:}...
    );
set(AX(2), 'YColor', [0.30,0.75,0.93]);
xlabel('$\mathrm{Time\ (sec)}$', 'Interpreter', 'latex');



%% Plot bifurcation diagram as function of G
nexttile;
lw = {'LineWidth', 2};
ltx={'Interpreter','latex'};
fna={'FontName','Courier New'};
fns={'FontSize',25};
% Set the remaining axes properties
set(gca,'FontSize',25,'LineWidth',...
    1);
% file location: nondim_equilibria_periodic orbits for G
nondimfile='bdeps0244';

nondimpo='POinput0';

coco_plot_bd(epplot, nondimfile, 'G','x',@(x)x(3,:));

grid on; box on; hold on;
% dim data
coco_plot_bd(thm, nondimpo, 'G','MAX(x)',@(x)x(3,:));
coco_plot_bd(thm, nondimpo, 'G','MIN(x)',@(x)x(3,:));

ylabel({'$y_3(exc)$'},'Interpreter','latex');

% Create xlabel
xlabel('G',fns{:},fna{:},ltx{:});
xlim([0 8]);

legend('$p=0$ stable', '$p=0$ unstable','', 'SNP','', '$p=120$', 'Interpreter','latex' );

ylabel('$\mathrm{Frequency\ (Hz)}$', 'Interpreter', 'latex');
%%
saveas(figure(1),'fig1.svg','svg');








