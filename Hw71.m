close all
clear all
clc
%Load the measurement data
S = load('Robotmes.mat');
Y = S.y;
t = 1:100;

%Initialize the system variables
A = [1,1,0;0,0,0;0,0,1];
U = [0;0.8;0];
Gamma = [0;0.1;0];
H = [1,0,0;-1,0,1];
T = diag([10 10]);
Q = 1;
X = [Y(1,1);0.8;Y(1,1)+Y(2,1)];
P = [10 0 10; 0 0.01 0; 10 0 20];

% Intialize the priors
Pos = [Y(1,1)];
Vel = [0.8];
Wall = [Y(1,1)+Y(2,1)];

% Additional vectors to stores individual variances
Pos_Var = [P(1,1)];
Vel_Var = [P(2,2)];
Wall_Var = [P(3,3)];

%Kalman Filter Loop

G=2;
while G<= length(Y)
%Prediction Step
X_pred = A*X + U;
P_pred = A*P*A' + Gamma*Q*Gamma';

%Calculate the kalman gain
K = P_pred*H'*inv(H*P_pred*H'+T);

%Update Step
X_upd = X_pred + K*(Y(:,G)-H*X_pred);
P_upd = (eye(3,3) - K*H)*P_pred;

%Store the estimated values
Pos = [Pos X_upd(1,1)];
Vel = [Vel X_upd(2,1)];
Wall = [Wall X_upd(3,1)];

%Store the variances
Pos_Var = [Pos_Var P_upd(1,1)];
Vel_Var = [Vel_Var P_upd(2,2)];
Wall_Var = [Wall_Var P_upd(3,3)];

%Updating the priors for the next step
X = X_upd;
P = P_upd; 
G = G+1;
end 


%reverse kalman 

% X_rev = X_upd;
% P_rev = P_upd;

X_rev = [Y(1,100);1.6;Y(1,100)+Y(2,100)];
P_rev = [10 0 10; 0 1 0; 10 0 20];
U_rev = [0;0;0];
Gamma_rev = [0;1;0];
% Intialize the priors
% Pos_rev = [Y(1,100)];
% Vel_rev = [0.8];
% Wall_rev = [Y(1,100)+Y(2,100)];
Pos_rev = [];
Vel_rev = [];
Wall_rev = [];

% Additional vectors to stores individual variances
Pos_Rev_Var = [P_rev(1,1)];
Vel_Rev_Var = [P_rev(2,2)];
Wall_Rev_Var = [P_rev(3,3)];


G=length(Y);

while G >= 1
%Prediction Step
X_rev_pred = A*X_rev + U_rev;
P_rev_pred = A*P_rev*A' + Gamma_rev*Q*Gamma_rev';

%Calculate the kalman gain
K_rev = P_rev_pred*H'*inv(H*P_rev_pred*H'+T);

%Update Step
X_rev_upd = X_rev_pred + K_rev*(Y(:,G)-H*X_rev_pred);
P_rev_upd = (eye(3,3) - K_rev*H)*P_rev_pred;

%Store the estimated values
Pos_rev = [Pos_rev X_rev_upd(1,1)];
Vel_rev = [Vel_rev X_rev_upd(2,1)];
Wall_rev =[Wall_rev X_rev_upd(3,1)];

%Store the variances
Pos_rev_Var = [Pos_Var P_upd(1,1)];
Vel_rev_Var = [Vel_Var P_upd(2,2)];
Wall_rev_Var = [Wall_Var P_upd(3,3)];

%Updating the priors for the next step
X_rev = X_rev_upd;
P_rev = P_rev_upd; 
G= G-1;
end 

figure;
plot(t,(Y(1,:)),'g-',t,Pos(1,:),'r-',t,fliplr(Pos_rev(1,:)),'b-','LineWidth',1.5);
ylabel('Magnitude of Position measurements and its Estimated value');
xlabel('Discrete time intervals k');
title('Plot of Position measurements and its Estimated values with time');
legend('Position measurement','Reverse Estimation','Forward Estimation','Location','southeast');

figure(2);
plot(t,Y(1,:)+Y(2,:),'g-',t,Wall(1,:),'r-',t,fliplr(Wall_rev(1,:)),'b-','LineWidth',1.5);
ylabel('Magnitude of Wall distance measureants and its Estimated Value');
xlabel('Discrete time intervals k');
title('Plot of Wall distance measureants and its Estimated Values with time');
legend('Wall distance measureant','Reverse Estimation','Forward Estimation','Location','southeast');

figure(3);
plot(t,Vel(1,:),'r-',t,fliplr(Vel_rev(1,:)),'b-','LineWidth',1.5);
ylabel('Magnitude of Velocity and its Variance');
xlabel('Discrete time intervals k');
title('Plot of Velocity and its variance with time');
legend('Forward Estimation', 'Reverse Estimation','Location','southeast');
