%% plot po_grazing_illu.pdf for thesis
clear
parnames={'jnew','G','P','alpha1','alpha2','alpha3','alpha4','b_star','c','beta',...
    'epsilon','x01','x02','x03','a1','a2','a3'};
format compact
cl=[parnames;num2cell(1:length(parnames))];
ip=struct(cl{:});                % automaticaly number parameters
newsigm=@(y,y0,epsilon)1./(1+exp(((y0)-y)./(epsilon)));
%% plot illustration
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
thm3=thm1;
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
%exportgraphics(fig,'../PhD-Thesis-draft/figures/po_grazing_illu.pdf','ContentType','vector');
