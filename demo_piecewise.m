clear
%%
% run startup of coco August 2023
%%
dimpar=load([pwd(),filesep,'dimeparam.mat']);% dowonload param scribt file 
c=dimpar.r*dimpar.v0;
beta=2*dimpar.a/(dimpar.B*dimpar.r*dimpar.C*2*dimpar.e0);% 0.024;;
alpha1=1;
alpha2=0.8;
alpha3=0.25;
alpha4=0.25;
y01=c*beta;
y02=c*beta/alpha3;                                                                                                                                                          
y03=c*beta/alpha1;
jnew=2; 
% read data to take value of G and 
fileread=coco_bd_read('po_smalleps001UZ');
labs=coco_bd_labs(fileread,'uz');
profile=po_read_solution('po_smalleps001UZ',labs(4)); % alpha
xnew=[profile.xbp(1:end-1,:);profile.xbp];
tnew=[profile.tbp(1:end-1);profile.tbp+profile.T];
G=profile.p(2); % G=1.7;
u1=[jnew/G; 0]; % input
b=0.5;
u2=[alpha4*jnew/b; 0]; % Example input
y3=alpha2*jnew/G;
y1thr=y3-y01;
%% time points I got them from fsolverfortimepoints.m
% load the result
t=load('initialdataoft');
T=t.t(1);
t_s1_on=t.t(2);
t_s1_off=t.t(3);
t_s2_on=t.t(4);
t_s2_off=t.t(5);

%%
% initial guess  y1_min and t1min read out from  data 'po_smalleps001UZ'
xnew=[profile.xbp(1:end-1,:);profile.xbp];
tnew=[profile.tbp(1:end-1);profile.tbp+profile.T];
plot(tnew, xnew(:,1));
%%
varnames={'y01','y02','y03','bstar','jnew','G','alpha2','alpha4',...
    'ts1off','ts1on','ts2off','ts2on','T',...
    't1min','y1_min','y1dot_min'};
format compact
cl=[varnames;num2cell(1:length(varnames))];
iv=struct(cl{:});                % automatically number parameters
u0=NaN(length(varnames),1);
u0( [iv.y01;iv.y02;iv.y03;iv.bstar;iv.jnew;iv.G; iv.alpha2;iv.alpha4])=...
    [ y01; y02;     y03;      b;          jnew;  G;  alpha2;  alpha4];
u0( [iv.ts1off;iv.ts1on;iv.ts2off;iv.ts2on;  iv.T])=...
    [t_s1_off;    t_s1_on   ;  t_s2_off  ; t_s2_on; T]; 
u0( [iv.t1min;iv.y1_min;iv.y1dot_min])=...
    [5.8       ;  0.137;     0       ]; % fix y1dot at zero
%% check pw_orbit
% create an empty prbolem 
prob=coco_prob();
[~,res]=pw_orbit(prob,iv,u0); % test my func(algebraic equations), res should be small 
prob=powbasic_coco(prob,iv,u0);
% add user point where 'y1_min'=y03
prob = coco_add_event(prob, 'graze', 'y1_min', u0(iv.y03));
%% get solutions of unknowns 
unkowns={'ts1on','ts2off','ts2on',...,
    'T','t1min','y1_min'};
bd_p=coco(prob,'initial_sol',[],0,unkowns); % this does the job as  fsolve 
%%
% we vary bstar and ymin along with time of switches  and fix all other paramter 
unkowns={'bstar','ts1on','ts2off',...,
         'ts2on','T','t1min',...,
          'y1_min'};
bd_p=coco(prob,'bstar_curve',[],1,unkowns,{[0.2,0.51]}); %  put 1 as this give us curve 

%%
% now take grazing bifurcatiom from 'bstar_curve'
% continue grazing in two param. b*,G 
upar=bd_getpoint('bstar_curve','graze',varnames);
prob=coco_prob();
[~,res]=pw_orbit(prob,iv,upar)
prob=powbasic_coco(prob,iv,upar);
%  vary bstar & G with fixed y1min
unkowns={'bstar','G', 'ts1on','ts2off','ts2on','T',...
    't1min'};
% prob = coco_set(prob, 'cont', 'PtMX',  [200 0 ],'NPR',10);
bd_p=coco(prob,'bstar_G',[],1,unkowns,{[0.2,0.51],[0.5 20]});

%% plotting curve of grazing bifurcation in b^*,G 

filename='bstar_G';
bda=coco_bd_read(filename);
u01=coco_bd_col(bda,varnames);
figure;hold on;
%%
plot(u01(iv.bstar,:),u01(iv.G,:),'Color','k')
plot(u01(iv.bstar,uzidx),u01(iv.y1_min,uzidx),'kd' , 'MarkerFaceColor' , 'g'  ,'MarkerSize' , 16);
xlabel('$bstar$');
ylabel('$G$');

