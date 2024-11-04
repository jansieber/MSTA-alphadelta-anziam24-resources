% run file symnondritjan_eps
%%
tic
clear;
parnames={'jnew','G','P','alpha1','alpha2','alpha3','alpha4','b_star','c','beta',...
    'epsilon','x01','x02','x03','a1','a2','a3'};
format compact
cl=[parnames;num2cell(1:length(parnames))];
ip=struct(cl{:});
F=sco_gen(@symnondritjan_eps);       % F and F2 are the same
funcs = {F(''),F('x'),F('p')};
dimpar=load('dimeparam.mat');
beta=2*dimpar.a/(dimpar.B*dimpar.r*dimpar.C*2*dimpar.e0); %0.024
alpha1=1; alpha3=0.25;
jnew=2;
c=dimpar.r*dimpar.v0;
% here call x01=y01 % x02=y02 % x03=y03
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
% Hopf data in bstar,G
ubg=bd_getpoint('GandbstarHB','UZ',parnames);
% Extract G and b_star values from data
G_values =ubg(ip.G,:);
% take  Uz 
b_star_values = ubg(ip.b_star,:);
dtG = 0.01: 0.1:2.5;
periodsarray=NaN(length(dtG),length(b_star_values));
G_mod=NaN(length(dtG),length(b_star_values));
% Loop over G_values and b_star_values
for j =1:length(b_star_values)
    G = G_values(j);
    for i = 1:length(dtG) 
        G_mod(i,j) = G + dtG(i); 
        % Set the parameters for the ODE
        par(ip.G) = G_mod(i,j);
        par(ip.b_star) = b_star_values(j);
        u0 = par'; 
        dim=7;
        x0=zeros(dim,1);
        % Solve the ODE with event detection
        f_ode = @(t, x) funcs{1}(x, u0);
        tspan = [0, 500];  
        options = odeset('Events', @(t, x) eventFunc(x, u0, ip), 'RelTol', 1e-6);
        sol = ode45(f_ode, tspan, x0, options);
        eventtime=(sol.xe(sol.xe > 300 & sol.ie == 1));% Find indices of events where ie == 1
        if length(eventtime) > 1
            % Calculate the period between successive events
            periods = diff(eventtime);  
            periodsarray(i,j)=mean(periods);
        end
    end
end
bstararry=repmat(b_star_values,length(dtG),1);
levels = [0, 1, 2, 5, 10, 20, 50, 100];
figure;
contourf(bstararry, G_mod, periodsarray,levels,'LineStyle','none');
save('bs_vals6050.mat', 'bstararry');
save('G_vals6050.mat', 'G_mod');
save('period_vals6050.mat','periodsarray');
 xlabel('$b^*$');
 ylabel('$G$');
toc
%%
function[check,terminal,direction]=eventFunc(y,u0,ip)
check=[y(4);y(1)-u0(ip.x03)]; %event condition
terminal=[0;0] ;
direction=[-1;-1];
end