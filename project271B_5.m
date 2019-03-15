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
    switchtime=0;stime=[];lambda=0.25;
    while switchtime<=9.99
        stime=[stime switchtime];
        switchtime=switchtime-1/lambda*log(1-unifrnd(0,1));
    end
    prob=unifrnd(0,1,1,length(stime)-1);
    
    if unifrnd(0,1)<=0.5
        at=100;
    else
        at=-100;
    end
    state0=[0;normrnd(0,200);at];
    state=state0;state_=[state0];
    s_estimate0=[0;0;0];
    s_estimate=s_estimate0;s_estimate_=[s_estimate0];
    error=[s_estimate0-state0];
    residual=[];
    index=2;temp=at;
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
        
        if index<=length(stime)
            if t<=stime(index)%&&prob(index-1)<=0.5
                state(3)=at;
                temp=state(3);
            else
                index=index+1;
                at=-at;
                state(3)=at;
                temp=state(3);
            end
      
        end
        d_state=F*state+G*wat;%dynamic propagate
        z=H*state+n;%measurenment
        state=state+d_state*dt;%state update
        state(3)=temp;
        state_=[state_ state];
        
        d_s_estimate=F*s_estimate+K*(z-H*s_estimate);%estimation update
        residual=[residual z-H*s_estimate];%calculate the residual
        s_estimate=s_estimate+d_s_estimate*dt;%state estimation update
        s_estimate_=[s_estimate_ s_estimate];
        error=[error s_estimate-state];%store the errors
    end

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




