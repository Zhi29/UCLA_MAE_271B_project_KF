Num=1000;
for i=1:Num
    project271B;%for problem 1-4
    project271B_5;%for problem 5
    if i==1 %define useful variables
    error_ave=zeros(3,length(error));
    error_store=zeros(3,length(error),Num);
    residual_store=zeros(1,length(residual),Num);
    var_sum=zeros(3,3);
    end
    error_ave=error_ave+error; % compute error sum
    error_store(:,:,i)=error;
    residual_store(:,:,i)=residual;
end
error_ave=error_ave/Num;  %compute e_ave

%% compute P_ave ======================
P_ave=[];
for tt=1:size(P_,2)/3
    for i=1:Num
        term=error_store(:,tt,i)-error_ave(:,tt);
        sum_p=term*term';
        var_sum=var_sum+sum_p;%3x3
    end
    var_sum=var_sum/(Num-1);
    P_ave=[P_ave var_sum];
    var_sum=zeros(3,3);
end
%% =======================================

%% residual orthogonality ==================
sum=0;
k=101;j=689;
        for i=1:Num
            sum=sum+residual_store(:,k,i)*residual_store(:,j,i)';
        end     
sum=sum/Num;
%% ==========================================

figure(5);
i=1:3:size(P_,2);
subplot 311
plot(error_ave(1,:),'linewidth',2);hold on;
plot(sqrt(P_(1,i)),'r--','linewidth',2);hold on;
plot(-sqrt(P_(1,i)),'r--','linewidth',2);grid on;
title('error average of position over realizations')
xlabel('time s');ylabel('e ave of position');
set(gca,'xticklabel',[10 9 8 7 6 5 4 3 2 1 0]);
xlim([0,1000]);
subplot 312
plot(error_ave(2,:),'linewidth',2);hold on;
plot(sqrt(P_(2,i+1)),'r--','linewidth',2);hold on;
plot(-sqrt(P_(2,i+1)),'r--','linewidth',2);grid on;
title('error average of velocity over realizations')
xlabel('time s');ylabel('e ave of velocity');
set(gca,'xticklabel',[10 9 8 7 6 5 4 3 2 1 0]);
xlim([0,1000]);
subplot 313
plot(error_ave(3,:),'linewidth',2);hold on;
plot(sqrt(P_(3,i+2)),'r--','linewidth',2);hold on;
plot(-sqrt(P_(3,i+2)),'r--','linewidth',2);grid on;
title('error average of acceleration over realizations')
xlabel('time s');ylabel('e ave of acceleration');
set(gca,'xticklabel',[10 9 8 7 6 5 4 3 2 1 0]);
xlim([0,1000]);

   
%% P variace and actual variance
figure(8);
ind=1:3:length(P_);
legend('Actual Error Variance from Monto Carlo Simulation','Priori Error Variance in Kalman Filter')
title('Variance Check for 1000 realizations for Pave and P');
subplot 331
plot(P_ave(1,ind),'k','linewidth',2);hold on; %%simulation and the average over an ensemble of realizations
plot(P_(1,ind),'r','linewidth',2); %% estimate variance of one realization
title('P(1,1)');
xlabel('time s');ylabel('error estimate variances');
set(gca,'xticklabel',[10 8 6 4 2 0]);
xlim([0,1000]);


subplot 332
plot(P_ave(1,ind+1),'k','linewidth',2);hold on; %%simulation and the average over an ensemble of realizations
plot(P_(1,ind+1),'r','linewidth',2); %% estimate variance of one realization
title('P(1,2)');
xlabel('time s');ylabel('error estimate variances');
set(gca,'xticklabel',[10 8 6 4 2 0]);
xlim([0,1000]);


subplot 333
plot(P_ave(1,ind+2),'k','linewidth',2);hold on; %%simulation and the average over an ensemble of realizations
plot(P_(1,ind+2),'r','linewidth',2); %% estimate variance of one realization
title('P(1,3)');
xlabel('time s');ylabel('error estimate variances');
set(gca,'xticklabel',[10 8 6 4 2 0]);
xlim([0,1000]);


subplot 334
plot(P_ave(2,ind),'k','linewidth',2);hold on; %%simulation and the average over an ensemble of realizations
plot(P_(2,ind),'r','linewidth',2); %% estimate variance of one realization
title('P(2,1)');
xlabel('time s');ylabel('error estimate variances');
set(gca,'xticklabel',[10 8 6 4 2 0]);
xlim([0,1000]);


subplot 335
plot(P_ave(2,ind+1),'k','linewidth',2);hold on; %%simulation and the average over an ensemble of realizations
plot(P_(2,ind+1),'r','linewidth',2); %% estimate variance of one realization
title('P(2,2)');
xlabel('time s');ylabel('error estimate variances');
set(gca,'xticklabel',[10 8 6 4 2 0]);
xlim([0,1000]);


subplot 336
plot(P_ave(2,ind+2),'k','linewidth',2);hold on; %%simulation and the average over an ensemble of realizations
plot(P_(2,ind+2),'r','linewidth',2); %% estimate variance of one realization
title('P(2,3)');
xlabel('time s');ylabel('error estimate variances');
set(gca,'xticklabel',[10 8 6 4 2 0]);
xlim([0,1000]);

subplot 337
plot(P_ave(3,ind),'k','linewidth',2);hold on; %%simulation and the average over an ensemble of realizations
plot(P_(3,ind),'r','linewidth',2); %% estimate variance of one realization
title('P(3,1)');
xlabel('time s');ylabel('error estimate variances');
set(gca,'xticklabel',[10 8 6 4 2 0]);
xlim([0,1000]);

subplot 338
plot(P_ave(3,ind+1),'k','linewidth',2);hold on; %%simulation and the average over an ensemble of realizations
plot(P_(3,ind+1),'r','linewidth',2); %% estimate variance of one realization
title('P(3,2)');
xlabel('time s');ylabel('error estimate variances');
set(gca,'xticklabel',[10 8 6 4 2 0]);
xlim([0,1000]);

subplot 339
plot(P_ave(3,ind+2),'k','linewidth',2);hold on; %%simulation and the average over an ensemble of realizations
plot(P_(3,ind+2),'r','linewidth',2); %% estimate variance of one realization
title('P(3,3)');
xlabel('time s');ylabel('error estimate variances');
set(gca,'xticklabel',[10 8 6 4 2 0]);
xlim([0,1000]);


%% P differrence variance plot
figure(6);
i=1:3:size(P_,2);
difference_var=P_ave-P_;
subplot 331
plot(difference_var(1,i));grid on;
title('P(1,1)');
xlabel('time s');ylabel('error estimate variances');
subplot 332
plot(difference_var(1,i+1));grid on;
title('P(1,2)');
xlabel('time s');ylabel('error estimate variances');
subplot 333
plot(difference_var(1,i+2));grid on;
title('P(1,3)');
xlabel('time s');ylabel('error estimate variances');
subplot 334
plot(difference_var(2,i));grid on;
title('P(2,1)');
xlabel('time s');ylabel('error estimate variances');
subplot 335
plot(difference_var(2,i+1));grid on;
title('P(2,2)');
xlabel('time s');ylabel('error estimate variances');
subplot 336
plot(difference_var(2,i+2));grid on;
title('P(2,3)');
xlabel('time s');ylabel('error estimate variances');
subplot 337
plot(difference_var(3,i));grid on;
title('P(3,1)');
xlabel('time s');ylabel('error estimate variances');
subplot 338
plot(difference_var(3,i+1));grid on;
title('P(3,2)');
xlabel('time s');ylabel('error estimate variances');
subplot 339
plot(difference_var(3,i+2));grid on;
title('P(3,3)');
xlabel('time s');ylabel('error estimate variances');
% 
% figure(7);
% plot(-100:100,sum2);
    
    
    
    
    

