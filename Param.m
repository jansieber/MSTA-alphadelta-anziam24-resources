%% model pramaters of dimensional version
clear; 
C=135; % connectivity constants
C1=C;
C2=0.8*C;
C3=0.25*C;
C4=0.25*C;
A=3.25;
B=22; % A,B mplitude of excitatory, inhibitory PSPs,resp.
a=100; b=50;% time constant of passive membrane 
P=0; % input to the model 
%% Param of sigmodial transformation
e0=2.5;% (2*e0=2*(2.5)=5)
r=0.56; %  the slpoe of sigma at v0(derivtive)
v0=6; %  value for which half(50%) of max  firing rate is attained                                                                                                                                                                                                                                                                                          
save('dimeparam.mat');