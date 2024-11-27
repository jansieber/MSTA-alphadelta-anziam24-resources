parnames={'jnew','G','P','alpha1','alpha2','alpha3','alpha4','b_star','c','beta',...
    'epsilon','x01','x02','x03','a1','a2','a3'};
format compact
cl=[parnames;num2cell(1:length(parnames))];
ip=struct(cl{:});
%a finer sampling over the range of bstar
bstargrid=linspace(0.2,0.5,500);
% HB in G,bstar
bG=bd_getpoint('GandbstarHB','UZ',parnames);
Ghopf1=interp1(bG(ip.b_star,:),bG(ip.G,:),bstargrid,'spline');
% sn in G,bstar
filename='SN_Gbstar_2_run';
ubg=bd_getpoint('SN_Gbstar_2_run','UZ',parnames);
Gsn1=interp1(ubg(ip.b_star,:),ubg(ip.G,:),bstargrid,'spline');

%%
save('Ghopf.mat', 'Ghopf');
save('Gsn.mat', 'Gsn');
%%
figure
plot(bstargrid,Ghopf1,'*r');
hold on
plot(bstargrid,Gsn1,'*b');
grid on;
