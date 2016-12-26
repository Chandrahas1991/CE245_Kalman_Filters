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

%Plots
plot(t,Y(1,:),'r-',t,Y(2,:),'g-','LineWidth',1.5);
ylabel('Magnitude of Position measurements and Wall measurements');
xlabel('Discrete time intervals k');
title('Plot of Position measurements and Wall measurements with time');
legend('Position measurement', 'Wall measurement','Location','northeast');

figure(2);
plot(t,Y(1,:),'g-',t,Pos(1,:),'r-','LineWidth',1.5);
ylabel('Magnitude of Position measurements and its Estimated value');
xlabel('Discrete time intervals k');
title('Plot of Position measurements and its Estimated value with time');
legend('Position measurement', 'Estimated value','Location','northeast');

figure(3);
plot(t,Y(1,:)+Y(2,:),'g-',t,Wall(1,:),'r-','LineWidth',1.5);
ylabel('Magnitude of Wall distance measureants and its Estimated Value');
xlabel('Discrete time intervals k');
title('Plot of Wall distance measureants and its Estimated Value with time');
legend('Wall distance measureant', 'Estimated Value','Location','northeast');

figure(4);
plot(t,Vel(1,:),'r-',t,Vel_Var(1,:),'g-','LineWidth',1.5);
ylabel('Magnitude of Velocity and its Variance');
xlabel('Discrete time intervals k');
title('Plot of Velocity and its variance with time');
legend('Velocity', 'Velocity variance','Location','northeast');

figure(5);
plot(t,Pos_Var(1,:),'r-',t,Wall_Var(1,:),'g-','LineWidth',1.5);
ylabel('Magnitude of Position Vairance and Wall distance variance');
xlabel('Discrete time intervals k');
title('Plot ofPosition Vairance and Wall distance variance with time');
legend('Position Vairance', 'Wall distance variance','Location','northeast');