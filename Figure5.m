% piece-wise periodic orbits  in the limit 
%% 
% figure 5 

%%
figure; nexttile

varnames={'y01','y02','y03','bstar','jnew','G','alpha2','alpha4',...
    'ts1off','ts1on','ts2off','ts2on','T',...
    't1min','y1_min','y1dot_min'};
format compact
cl=[varnames;num2cell(1:length(varnames))];
iv=struct(cl{:});      
%%

filename='bstar_curve';

bG=coco_bd_read(filename);   
u0=coco_bd_col(bG,varnames);
% G=1.7, bstar=0.5
i=41;

b=[1;u0(iv.bstar,i)];
u=[0;u0(iv.jnew,i)/u0(iv.G,i)];
v=[0;u0(iv.alpha4)*u0(iv.jnew)/u0(iv.bstar,i)];

%%
% use the formula 
y1per=yper(u0(iv.ts1on,i)-u0(iv.ts1off,i),u0(iv.T,i),b(1),u(2));
y2per=yper(u0(iv.ts2on,i)-u0(iv.ts2off,i),u0(iv.T,i),b(2),v(2));
tcy1=linspace(u0(iv.ts1off,i),u0(iv.ts1on,i),10000);
tcy1interval2=linspace(u0(iv.ts1on,i),u0(iv.ts1off,i)+u0(iv.T,i), 10000);
tcy2=linspace(u0(iv.ts2off,i),u0(iv.ts2on,i), 10000);
tcy2interval2=linspace(u0(iv.ts2on,i),(u0(iv.ts2off,i)+u0(iv.T,i)), 10000);
%%
% y1 component 

y1=pw_advance(y1per,tcy1-u0(iv.ts1off,i),b(1),u(1));
y1interval2=pw_advance(y1(:,end),tcy1interval2-u0(iv.ts1on,i),b(1),u(2));

%%

% y2 component 
y2=pw_advance(y2per,tcy2-u0(iv.ts2off,i),b(2),v(1));
y2interval2=pw_advance(y2(:,end),tcy2interval2-u0(iv.ts2on,i),b(2),v(2));

%%
hold on;
lwidth=3.5; 
plot(tcy1,y1(1,:),'LineWidth',lwidth);
plot(tcy1interval2,y1interval2(1,:),'LineWidth',lwidth,'Color',[0.64,0.08,0.18]);

plot(tcy2,y2(1,:),'LineWidth',lwidth);
plot(tcy2interval2,y2interval2(1,:),'LineWidth',lwidth);
grid on;
u0(iv.y1_min,i)
xline(u0(iv.t1min,i))
yline(u0(iv.y03,i))
plot(u0(iv.t1min,i),u0(iv.y1_min,i),'kd' , 'MarkerFaceColor' , 'g'  ,'MarkerSize' , 10);
%%

nexttile
% panel b, ynim=y03
% G=1.7, bstar=0.44
u0=bd_getpoint('bstar_curve','graze',varnames);
i=1;
b=[1;u0(iv.bstar,i)];
u=[0;u0(iv.jnew,i)/u0(iv.G,i)];
v=[0;u0(iv.alpha4)*u0(iv.jnew)/u0(iv.bstar,i)];
%%

% use the formula 
y1per=yper(u0(iv.ts1on,i)-u0(iv.ts1off,i),u0(iv.T,i),b(1),u(2));
y2per=yper(u0(iv.ts2on,i)-u0(iv.ts2off,i),u0(iv.T,i),b(2),v(2));
tcy1=linspace(u0(iv.ts1off,i),u0(iv.ts1on,i),10000);
tcy1interval2=linspace(u0(iv.ts1on,i),u0(iv.ts1off,i)+u0(iv.T,i), 10000);
tcy2=linspace(u0(iv.ts2off,i),u0(iv.ts2on,i), 10000);
tcy2interval2=linspace(u0(iv.ts2on,i),(u0(iv.ts2off,i)+u0(iv.T,i)), 10000);
%%
% y1 component 

y1=pw_advance(y1per,tcy1-u0(iv.ts1off,i),b(1),u(1));
y1interval2=pw_advance(y1(:,end),tcy1interval2-u0(iv.ts1on,i),b(1),u(2));

%%

% y2 component 
y2=pw_advance(y2per,tcy2-u0(iv.ts2off,i),b(2),v(1));
y2interval2=pw_advance(y2(:,end),tcy2interval2-u0(iv.ts2on,i),b(2),v(2));

%%
hold on;
lwidth=3.5; 
plot(tcy1,y1(1,:),'LineWidth',lwidth);
plot(tcy1interval2,y1interval2(1,:),'LineWidth',lwidth,'Color',[0.64,0.08,0.18]);

plot(tcy2,y2(1,:),'LineWidth',lwidth);
plot(tcy2interval2,y2interval2(1,:),'LineWidth',lwidth);

grid on;
u0(iv.y1_min,i)
xline(u0(iv.t1min,i))
yline(u0(iv.y03,i))
plot(u0(iv.t1min,i),u0(iv.y1_min,i),'kd' , 'MarkerFaceColor' , 'g'  ,'MarkerSize' , 10);

%%
function y=yper(toff2on,T,b,u)
M=@(t)pw_lin(t,b);
Id=eye(2);
y=(Id-M(T))\(Id-M(T-toff2on))*[u;0];
end
