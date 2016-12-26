close all
clear all
clc
%Load the measurement data
S = load('Robotmes.mat');
Y = S.y;
t = 1:100;

%Initialize the system variables
Phi = [1,1,0;0,1,0;0,0,1];
U = [0;0;0];
Gamma = [0;0.1;0];
H = [1,0,0;-1,0,1];
T = diag([10 10]);
Q = 1;
X = [Y(1,1);0;Y(1,1)+Y(2,1)];
P = [10 0 10; 0 1 0; 10 0 20];
X_pred_hist = [];
P_pred_hist = [];
X_upd_hist = [];
P_upd_hist =[];
% Intialize the priors
Pos = [];
Vel = [];
Wall = [];

% Additional vectors to stores individual variances
Pos_Var_Kal = [];
Vel_Var_Kal = [];
Wall_Var_Kal = [];

%Kalman Filter Loop

for G = 1:100
%Prediction Step
X_pred = Phi*X + U;
P_pred = Phi*P*Phi' + Gamma*Q*Gamma';
X_pred_hist = [X_pred_hist X_pred];
P_pred_hist = [P_pred_hist P_pred];

%Calculate the kalman gain
K = P_pred*H'*inv(H*P_pred*H'+T);

%Update Step
X_upd = X_pred + K*(Y(:,G)-H*X_pred);
P_upd = (eye(3,3) - K*H)*P_pred;
X_upd_hist = [X_upd_hist X_upd];
P_upd_hist = [P_upd_hist P_upd];


%Store the estimated values
Pos =  [Pos X_upd(1,1)];
Vel =  [Vel X_upd(2,1)];
Wall = [Wall X_upd(3,1)];

%Store the variances
Pos_Var_Kal =  [Pos_Var_Kal P_upd(1,1)];
Vel_Var_Kal =  [Vel_Var_Kal P_upd(2,2)];
Wall_Var_Kal = [Wall_Var_Kal P_upd(3,3)];

%Updating the priors for the next step
X = X_upd;
P = P_upd; 
G = G+1;
end
Kal_Final_Var = P_upd
% RTS Smoother 

X_est = X_upd;
P_est = P_upd;
X_est_hist = [X_est];
P_est_hist = [P_est];
M = length(Y)-1;
N = M+1;
Pos_Var_RTS =[P_est(1,1)];
Vel_Var_RTS =[P_est(2,2)];
Wall_Var_RTS =[P_est(3,3)];
while M>=1
    A = P_upd_hist(:,3*M-2:3*M)*Phi'*inv(P_pred_hist(:,(3*N)-2:3*N));
    X_est = X_upd_hist(:,M)+ A*(X_est - X_pred_hist(:,N));
    P_est = P_upd_hist(:,3*M-2:3*M) + A*(P_est -P_pred_hist(:,3*N-2:3*N))*A';
    X_est_hist = [X_est_hist X_est];
    P_est_hist = [P_est_hist P_est];
    Pos_Var_RTS = [Pos_Var_RTS P_est(1,1)];
    Vel_Var_RTS = [Vel_Var_RTS P_est(2,2)];
    Wall_Var_RTS =[Wall_Var_RTS P_est(3,3)];
    M = M-1;
    N = N-1;    
end

RTS_Final_Var = P_est


%Plots

figure(2);
plot(t,Y(1,:),'g-',t,Pos(1,:),'r-',t,fliplr(X_est_hist(1,:)),'b-','LineWidth',1.5);
ylabel('Magnitude of Position measurements and its Estimated value');
xlabel('Discrete time intervals k');
title('Plot of Position measurements and its Estimated value with time');
legend('Position Measurements','Forward Kalman','RTS Estimate','Location','southeast');

figure(3);
plot(t,Y(1,:)+Y(2,:),'g-',t,Wall(1,:),'r-',t,fliplr(X_est_hist(3,:)),'b-','LineWidth',1.5);
ylabel('Wall distance measureants and its Estimated Value');
xlabel('Discrete time intervals k');
title('Plot of Wall distance measureants and its Estimated Value with time');
legend('Wall distance measureant','Forward Kalman','RTS Estimate','Location','northeast');

figure(4);
plot(t,Vel(1,:),'r-',t,fliplr(X_est_hist(2,:)),'b-','LineWidth',1.5);
ylabel('Magnitude of Velocity ');
xlabel('Discrete time intervals k');
title('Plot of Velocity for Kalman Filter and RTS Smoother');
legend('Forward Kalman', 'RTS Estimate','Location','northeast');

figure(5);
plot(t,Pos_Var_Kal(1,:),'r-',t,Pos_Var_RTS(1,:),'b-','LineWidth',1.5);
ylabel('Magnitude of Position Variance');
xlabel('Discrete time intervals k');
title('Plot of Position Vairance for Kalman Filter and RTS Smoother');
legend('Kalman Position Variance', 'RTS Position Variance','Location','northeast');

figure(6);
plot(t,Wall_Var_Kal(1,:),'r-',t,Wall_Var_RTS(1,:),'b-','LineWidth',1.5);
ylabel('Magnitude of Wall measurement Variance');
xlabel('Discrete time intervals k');
title('Plot of Wall measurement Variance for Kalman Filter and RTS Smoother');
legend('Kalman Position Variance', 'RTS Position Variance','Location','northeast');

figure(7);
plot(t,Vel_Var_Kal(1,:),'r-',t,Vel_Var_RTS(1,:),'b-','LineWidth',1.5);
ylabel('Magnitude of Velocity Variance');
xlabel('Discrete time intervals k');
title('Plot of Velocity Variance for Kalman Filter and RTS Smoother');
legend('Kalman Velocity Variance', 'RTS Velocity Variance','Location','northeast');

