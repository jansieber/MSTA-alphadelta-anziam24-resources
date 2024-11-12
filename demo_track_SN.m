varnames={'bstar','jnew','G','alpha2','alpha4',...
    'a1','a2','a3','y01','y02','y03','y1','epsilon'...
    'eqb','eqbd'};
format compact
cl=[varnames;num2cell(1:length(varnames))];
iv=struct(cl{:});
% track saddle node of xeq
%choose pramter for grazing
u0=NaN(length(varnames),1);
u0( [iv.bstar;iv.jnew;  iv.G; iv.alpha2;iv.alpha4])=...
    [ 0.5;     2;        1.7;      0.8;   0.25  ];
u0( [iv.a1;iv.a2;iv.a3])=...
    [1;    0.25   ;  1];
u0([iv.epsilon;iv.y01;iv.y02;iv.y03])=...
    [0.001     ;0.084  ;0.336; 0.084];
u0([iv.y1])=0.001;
u0([iv.eqb])=0;
u0([iv.eqbd])=0;
uxeq=u0;
prob=coco_prob();
[~,res]=singularxeq(prob,iv,uxeq)
% update the value of eqb, eqbd
uxeq([iv.eqb,iv.eqbd])=res;
[~,res]=singularxeq(prob,iv,uxeq)
% call snbasic to create problem
% pass equation and its paramters to coco
prob=snbasic_coco(iv,uxeq);
% add user point
prob = coco_add_event(prob, 'inig', 'eqb', 0);
%% 
unknowns={'y1','eqb','eqbd','y01',  'y03', 'G','y02','bstar'};
coco(prob,'ini_y01Gsmall',[],1,unknowns,{[-1,1],[-1,1],[-20,20]});
%%
nodeini=bd_getpoint('ini_y01Gsmall','inig',varnames);
prob=snbasic_coco(iv,nodeini);
prob = coco_add_event(prob, 'inisn', 'eqbd', 0);
unknowns={'y1','eqbd','y01','y03',  'G','y02','bstar'};
coco(prob,'ini_eqGsmall',[],1,unknowns,{[-1,1],[-20,20]});

%%
%  track sn(y1) in eps,y01,y03
sn_ini=bd_getpoint('ini_eqGsmall','inisn',varnames);
prob=snbasic_coco(iv,sn_ini);
unknowns={'y1','epsilon','y01','y03',    'y02','bstar'};
% fix y03=0.084
prob = coco_add_event(prob, 'fixy03', 'y03', 0.084);
coco(prob,'sn_epsandy03',[],1,unknowns,{[-1,2],[0.0001,0.05],[0.0001,0.09]});
%%
%  slice of differnt epsilon  with fixed y03=0.084
sn_ini=bd_getpoint('sn_epsandy03','fixy03',varnames);
prob=snbasic_coco(iv,sn_ini);
unknowns={'y1','G','epsilon','y01',      'y03','y02','bstar'};
%
prob = coco_add_event(prob, 'diffepsvalue', 'epsilon',[0.01,0.005,0.015,0.02])
coco(prob,'snwithdiffeps2',[],1,unknowns,{[-1,2],[0.0001,20],[0.001,0.025]});
%%
sn_ini=bd_getpoint('snwithdiffeps2','diffepsvalue',varnames);
prob=snbasic_coco(iv,sn_ini(:,2));%sn_ini(:,1)
% free paramters, y01 is inactive as we gule y01 & y03
unknowns={'y1','G','bstar','y01',   'y03','epsilon','y02'};
coco(prob,'sn_Gbstareps015',[],1,unknowns,{[-1,2],[0.003,3],[0.1,0.5]});
%coco(prob,'sn_Gbstareps01',[],1,unknowns,{[-1,2],[0.003,3],[0.1,0.5]});
