figure; clf; 
nexttile;
epplot=ep_plot_theme('ep');
epplot = struct('special', {{'SN', 'HB'}});
epplot.SN={'kd' , 'MarkerFaceColor' , 'g'  ,'MarkerSize' , 14};
epplot.HB={'kd' , 'MarkerFaceColor' , 'r'  ,'MarkerSize' , 14};
epplot.lspec{1}={'b-' ,   'LineWidth',   2};
epplot.lspec{2}={'b--'  ,  'LineWidth' ,   2};
coco_plot_bd(epplot,'smalleps=001' , 'G','x',@(x)x(1,:));
thm=po_plot_theme('po');
thm = struct('special', {{'SN'}});
thm.SN={'kd' , 'MarkerFaceColor' , 'g'  ,'MarkerSize' , 14};
thm.lspec{1}={'k-' ,   'LineWidth',   2};
thm.lspec{2}={'k--'  ,  'LineWidth' ,   2};
hold on; 
coco_plot_bd(thm, 'po_smalleps001UZ', 'G','MAX(x)',@(x)x(1,:))
coco_plot_bd(thm, 'po_smalleps001UZ', 'G','MIN(x)',@(x)x(1,:))
ylabel('$\mathrm{y_1}$(pc)', 'Interpreter', 'latex');
xlabel('$\mathrm{G}$', 'Interpreter', 'latex');        
%%
% y2 
nexttile
coco_plot_bd(epplot,'smalleps=001' , 'G','x',@(x)x(2,:));
hold on; 
coco_plot_bd(thm, 'po_smalleps001UZ', 'G','MAX(x)',@(x)x(2,:))
coco_plot_bd(thm, 'po_smalleps001UZ', 'G','MIN(x)',@(x)x(2,:))
ylabel('$\mathrm{y_2}$(inh)', 'Interpreter', 'latex');
xlabel('$\mathrm{G}$', 'Interpreter', 'latex');      
%%
% y3
coco_plot_bd(epplot,'smalleps=001' , 'G','x',@(x)x(3,:));
hold one; 
coco_plot_bd(thm, 'po_smalleps001UZ', 'G','MAX(x)',@(x)x(3,:))
coco_plot_bd(thm, 'po_smalleps001UZ', 'G','MIN(x)',@(x)x(3,:))
ylabel(ax,'$y_3$(exc)',ltx{:});
xlabel(ax,'$G$',ltx{:});
xlim([1.3,2.3]);
