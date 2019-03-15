close all;
clear all;
clc;
rng('shuffle');


%% initials
tau=2;%sec
Vc=300;%300 ft/sec
tf=10;%sec
R1=15e-6;%rad^2.sec
R2=1.67e-3;%rad^2.sec^3

F=[0,1,0;0,0,-1;0,0,-1/tau];
B=[0;1;0];
G=[0;0;1];

W=[0,0,0;0,0,0;0,0,10000];

P0=[0,0,0;0,40000,0;0,0,10000];
P=P0;
K_=[];P_=P0;dt=0.01;

%% dynamic simulation and Kalman Filter
state0=[0;normrnd(0,200);normrnd(0,100)];
state=state0;state_=[state0];
s_estimate0=[0;0;0];
s_estimate=s_estimate0;s_estimate_=[s_estimate0];
error=[s_estimate0-state0];
residual=[];

for t=0:0.01:9.99
    H=[1/(Vc*(tf-t)),0,0];
    V=R1+R2/(tf-t)^2;
    
    d_P=F*P+P*F'-1/V*P*H'*H*P+W; %Riccati Equation
    K=P*H'*inv(V); %Kalman Gain
    K_=[K_ K];
    P=d_P*dt+P; % Propagate P
    P_=[P_ P];
    
    wat=normrnd(0,sqrt(10000/dt));%process noise 
    n=normrnd(0,sqrt(V/dt));%measurement noise
    d_state=F*state+G*wat;%dynamic propagate
    z=H*state+n;%measurenment
    state=state+d_state*dt;%state update
    state_=[state_ state];
        
    d_s_estimate=F*s_estimate+K*(z-H*s_estimate);%estimation update
    residual=[residual z-H*s_estimate];%calculate the residual
    s_estimate=s_estimate+d_s_estimate*dt;%state estimation update
    s_estimate_=[s_estimate_ s_estimate];
    error=[error s_estimate-state];%store the errors
end

figure(1);
plot(K_(1,:),'b','Linewidth',2);hold on
plot(K_(2,:),'r-.','Linewidth',2);hold on;
plot(K_(3,:),'k--','Linewidth',2);grid on;
title('Filter Gain history for missile intercept estimation');
legend('K1','K2','K3');
xlabel('time-to-go/sec');ylabel('Kalman Filter Gain');
set(gca,'xticklabel',[10 9 8 7 6 5 4 3 2 1 0]);

figure(2);
i=1:3:size(P_,2);
plot(sqrt(P_(1,i)),'b','Linewidth',2);hold on;
plot(sqrt(P_(2,i+1)),'r-.','Linewidth',2);hold on;
plot(sqrt(P_(3,i+2)),'k--','Linewidth',2);grid on;
title('Evolution of the Estimation Error RMS for Missile Intercept Example');
legend('RMS error in position','RMS error in velocity','RMS error in acceleration');
xlabel('time-to-go/sec');ylabel('Standard deviation of the state error');
set(gca,'xticklabel',[10 9 8 7 6 5 4 3 2 1 0]);
xlim([0,1000]);

figure(3);
subplot 311
plot(state_(1,:),'b','Linewidth',2);hold on;
plot(s_estimate_(1,:),'r--','Linewidth',2);grid on;
title('relative vertical postion of pursuer with respect to the target');
legend('state','state estimation');
xlabel('time-to-go/sec');ylabel('position ft');
set(gca,'xticklabel',[10 9 8 7 6 5 4 3 2 1 0]);
xlim([0,1000]);
subplot 312
plot(state_(2,:),'b','Linewidth',2);hold on;
plot(s_estimate_(2,:),'r--','Linewidth',2);grid on;
title('relative vertical velocity of pursuer with respect to the target');
legend('state','state estimation');
xlabel('time-to-go/sec');ylabel('velocity ft/sec');
set(gca,'xticklabel',[10 9 8 7 6 5 4 3 2 1 0]);
xlim([0,1000]);
subplot 313
plot(state_(3,:),'b','Linewidth',2);hold on;
plot(s_estimate_(3,:),'r--','Linewidth',2);grid on;
title('relative vertical acceleration of pursuer with respect to the target');
legend('state','state estimation');
xlabel('time-to-go/sec');ylabel('acceleration ft/sec^2');
set(gca,'xticklabel',[10 9 8 7 6 5 4 3 2 1 0]);
xlim([0,1000]);

figure(4);
subplot 311
plot(error(1,:),'LineWidth',2);hold on;
plot(sqrt(P_(1,i)),'r--','LineWidth',2);hold on;
plot(-sqrt(P_(1,i)),'r--','LineWidth',2);grid on;
xlabel('time-to-go/sec');ylabel('position/ft');
title('error in the Kalman Filter estimates for single realization');
set(gca,'xticklabel',[10 9 8 7 6 5 4 3 2 1 0]);
xlim([0,1000]);
subplot 312
plot(error(2,:),'LineWidth',2);hold on;
plot(sqrt(P_(2,i+1)),'r--','LineWidth',2);hold on;
plot(-sqrt(P_(2,i+1)),'r--','LineWidth',2);grid on;
xlabel('time-to-go/sec');ylabel('position/ft');
title('error in the Kalman Filter estimates for single realization');
set(gca,'xticklabel',[10 9 8 7 6 5 4 3 2 1 0]);
xlim([0,1000]);
subplot 313
plot(error(3,:),'LineWidth',2);hold on;
plot(sqrt(P_(3,i+2)),'r--','LineWidth',2);hold on;
plot(-sqrt(P_(3,i+2)),'r--','LineWidth',2);grid on;
xlabel('time-to-go/sec');ylabel('position/ft');
title('error in the Kalman Filter estimates for single realization');
set(gca,'xticklabel',[10 9 8 7 6 5 4 3 2 1 0]);
xlim([0,1000]);


