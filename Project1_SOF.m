close all
clear all
clc

XY = csvread('XYData_cm.csv');
HA = csvread('HeadingAngle_rad.csv');
dT = 0.33;
V0 = sqrt((XY(2,1)-XY(1,1))^2 + (XY(2,2)-XY(1,2))^2)/dT;
X = [XY(1,1); XY(1,2);-V0;pi];
sX = 0.2;
sY = 0.2;
gamma = [0,0;0,0;sX*sqrt(dT),0;0,sY*sqrt(dT)];
Q = [1,0;0,1];
R = [0.0588,0;0,0.0588];
H = [1,0,0,0;0,1,0,0];
X_est = [X(1,1)];
Y_est = [X(2,1)];
vel_hist = [0];
Vel_est =[0];
Theta_est =[X(4,1)];
P = [sX,0,0,0;
     0,sY,0,0;
     0,0,1.07,0;
     0,0,0,0.0304];
 P11 =[P(1,1)];
 P22 =[P(2,2)];
 P33 =[P(3,3)];
 P44 =[P(4,4)];
 P11_hist =[P11];
 P22_hist =[P22];
 P33_hist =[P33];
 P44_hist =[P44];
 b =[];
 b1=[];
 b2=[];

for t = 2:length(XY)    
% Calulating the Bias terms for the prediction step. 
%         jacobian([1;0;dT*cos(x4);-dT*x3*sin(x4)],[x1 x2 x3 x4])
           b1 = [ 0, 0, 0,  0;
                0, 0, 0,  0;
                0, 0, 0, -dT*sin(X(4,1));
                0, 0, -dT*sin(X(4,1)), -dT*X(3,1)*cos(X(4,1))];
%             jacobian([0;1;dT*sin(x4);dT*x3*cos(x4)],[x1 x2 x3 x4])
            b2 = [ 0, 0, 0,  0;
                   0, 0, 0,  0;
                   0, 0, 0, dT*cos(X(4,1));
                   0, 0, dT*cos(X(4,1)), -dT*X(3,1)*sin(X(4,1))];
               
            b = 0.5*[trace(b1*P);trace(b2*P);0;0];     
    
    %Prediction Step
    x1 = X(1,1) + X(3,1)*dT*cos(X(4,1));
    x2 = X(2,1) + X(3,1)*dT*sin(X(4,1));
    x3 = X(3,1);
    x4 = wrapToPi(X(4,1));
    
    X_Pred = [x1;x2;x3;x4] + b;
    
%   jacobian([x1+(x3*dT*cos(x4));x2+(dT*x3*sin(x4));x3;x4],[x1 x2 x3 x4])
    Phi = [ 1, 0, dT*cos(x4), -x3*dT*sin(x4);
            0, 1, dT*sin(x4),  x3*dT*cos(x4);
            0, 0,       1,           0;
            0, 0,       0,           1];
        
    P_Pred = Phi*P*Phi' + gamma*Q*gamma';
    
    
        
    
    P_upd = P_Pred - P_Pred*H'*inv(H*P_Pred*H'+R)*H*P_Pred;
    K = P_upd*H'*inv(R);
    X_upd = X_Pred + K*(XY(t,:)'-H*X_Pred);
    X_upd(4,1) = wrapToPi(X_upd(4,1));
    
    X_est = [X_est X_upd(1,1)];
    Y_est = [Y_est X_upd(2,1)];
    Vel_est = [Vel_est X_upd(3,1)];
    Theta_est = [Theta_est X_upd(4,1)];
    P11_hist =[P11_hist P_upd(1,1)];
    P22_hist =[P22_hist P_upd(2,2)];
    P33_hist =[P33_hist P_upd(3,3)];
    P44_hist =[P44_hist P_upd(4,4)];
    
    vel = sqrt((XY(t,1)-XY(t-1,1))^2 + (XY(t,2) - XY(t-1,2))^2)/dT;
    
    vel_hist = [vel_hist vel ];
    
    X = X_upd;
    P = P_upd; 
    
end


t= 1:length(XY);
%Plots
figure(1);
plot(t,XY(:,1),'r-',t,X_est(1,:),'g-');
ylabel('Magnitude of X Position measurements in cm');
xlabel('Discrete time intervals k in seconds');
title('Plot of X Position measurements and Estimates with time');
legend('measurement', 'Estimate','Location','northeast');

figure(2);
plot(t,XY(:,2),'r-',t,Y_est(1,:),'g-');
ylabel('Magnitude of Y Position measurements in cm');
xlabel('Discrete time intervals k in seconds');
title('Plot of Y Position measurements and Estimates with time');
legend('measurement', 'Estimated value','Location','northeast');

figure(3);
plot(XY(:,1),XY(:,2),'r-','LineWidth',1.5);
hold on;
plot(X_est,Y_est,'g-','LineWidth',1.5);
ylabel('Y Position measurements in cm');
xlabel('X Position measurements in cm');
title('Plot of Position measurements and Estimate with time');
legend('measurement', 'Estimate','Location','northeast');

figure(4);
plot(t,HA(:,1),'r-',t,Theta_est(1,:),'g-','LineWidth',1.5);
ylabel('Magnitude of Heading Angle in radians');
xlabel('Discrete time intervals k in seconds');
title('Plot of Heading Angle and its estimates with time');
legend('measurement', 'Estimate','Location','northeast');


figure(5);
plot(t,vel_hist,'r-',t,Vel_est(1,:),'g-','LineWidth',1.5);
ylabel('Magnitude of X Position measurements in cm/sec');
xlabel('Discrete time intervals k in seconds');
title('Plot of Velocity Estimate with time');
legend('measurement','Estimate', 'Location','southeast'); 

figure(6);
plot(t,P11_hist,t,P22_hist,t,P33_hist,t,P44_hist);
ylabel('Magnitude of variance');
xlabel('Discrete time intervals k in seconds');
title('Plot of State Variance with time');
legend('P11','P22','P33','P44', 'Location','northeast');
    
    
    