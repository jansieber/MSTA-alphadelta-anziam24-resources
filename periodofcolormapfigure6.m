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
beta=2*dimpar.a/(dimpar.B*dimpar.r*dimpar.C*2*dimpar.e0);% 0.024;; 
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
 load('Ghopf.mat');
 load('Gsn.mat');
dt=0.01;
m=60;

bstargrid=linspace(0.2,0.5,50);
nbs=length(bstargrid);
bs_vals=repmat(bstargrid,m,1); % assuming bstargrid is 1 x nbs

% Preallocate arrays for storing the results (Periods)
period_vals=NaN(m,nbs);

G_vals=NaN(m,nbs);
for j=1:nbs
 G_vals(:,j)=linspace(Ghopf(j)+dt,Gsn(j)-dt,m);
end
%%
 % j for column, i for rows 
 % fix j and go over  rows 

for j=1:nbs
   
    for i =1:m
        par(ip.G) = G_vals(i,j);
        par(ip.b_star) = bs_vals(i,j); % Keep b_star fixed
        % Take the corresponding b_star value
        u0 = par';% Update parameters
        dim=7;
        x0=zeros(dim,1);
        % Solve the ODE with event detection
        f_ode = @(t, x) funcs{1}(x, u0);
        tspan = [0, 500];  % Make sure the time span is long enough to capture multiple events
        options = odeset('Events', @(t, x) eventFunc(x, u0, ip), 'RelTol', 1e-6);
        sol = ode45(f_ode, tspan, x0, options);
        % create logical array [ture,false]
        % Extract event times for the first event type (or modify for other events)
        eventtime=(sol.xe(sol.xe > 300 & sol.ie == 1));% Find indices of events where ie == 1
        if length(eventtime) > 1
            % Calculate the period between successive events
            periods = diff(eventtime);  % Calculate the differences (i.e., periods)
            %disp('Periods of events:');
            %disp(periods);
            %disp([i,j]);
            period_vals(i,j)=mean(periods);
  if j == 40
        fprintf('j=%d\n', j);
    end

    end
    end
end

% use contour levels  to distinguish  alpha and delta

 figure;
   contourf(bs_vals, G_vals, period_vals,'LineStyle','none');
  save('bs_vals6050.mat', 'bs_vals');
  save('G_vals6050.mat', 'G_vals');
  save('period_vals6050.mat', 'period_vals');
 toc

%%
% takes single value y(7,1) at each time points
function[check,terminal,direction]=eventFunc(y,u0,ip)
check=[y(4);y(1)-u0(ip.x03)]; %event condition
terminal=[0;0] ;% does not stop when event happen
direction=[-1;-1];% derivative change from + to -
end

