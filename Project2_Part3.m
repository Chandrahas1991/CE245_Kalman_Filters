close all
clear all
clc

D1 = load('dataset3.mat');

%Dataset 1 Data
dt = 0.1;
pitch = degtorad(D1.pitch);
roll = degtorad(D1.roll);
yaw = degtorad(D1.yaw);
xm = D1.xm;
ym = D1.ym;
vx(1) = (xm(2)-xm(1))/dt(1);
vy(1) = (ym(2)-ym(1))/dt(1);
Y = [xm;ym;yaw];
thrust = D1.thrust;

g =9.81;
sX = 0.01;
sY = 0.01;
sXd = 2*sX/(dt(1)^2);
sYd = 2*sY/(dt(1)^2);
sYaw = degtorad(0.01);

P =[];
Q = [1,0,0;...
     0,1,0;...
     0,0,1];
R = [0.01^2,0,0;...
     0,0.01^2,0;...
     0,0,(pi/180)^2];
H = [1,0,0,0,0;...
     0,1,0,0,0;...
     0,0,0,0,1];
X  = [xm(1);ym(1);vx(1);vy(1);yaw(1)];
Xm = [xm(1)];
Ym = [ym(1)];
Yaw_est = [yaw(1)];

P = [sX,0,0,0,0;
     0,sY,0,0,0;
     0,0,sXd,0,0;
     0,0,0,sYd,0;
     0,0,0,0,sYaw];
 
 
 for i = 2:length(xm)
    X(1) = X(1)+ dt*X(3);
    X(2) = X(2)+ dt*X(4);
    X(3) = X(3)+ dt*g*(sin(X(5))*sin(roll(i))+cos(X(5))*cos(roll(i))*sin(pitch(i)));
    X(4) = X(4)+ dt*g*(sin(X(5))*cos(roll(i))*sin(pitch(i))-cos(X(5))*sin(roll(i)));
    X(5) = X(5);
    
    X_Pred = X;
%   jacobian([X1+ dtrec(i)*X4,X2+ dtrec(i)*X5,X3+ dtrec(i)*X6,X4+ dtrec(i)*g*(sin(X7)*sin(roll(i))+cos(X7)*cos(roll(i))*sin(pitch(i))),X5+ dtrec(i)*g*(sin(X7)*cos(roll(i))*sin(pitch(i))-cos(X7)*sin(roll(i))),X6+ dtrec(i)*g*cos(roll(i))*cos(pitch(i))-g, X7],[X1 X2 X3 X4 X5 X6 X7])
 
    Phi =[ 1, 0, dt,  0, 0;...
           0, 1,  0, dt, 0;...
           0, 0,  1,  0, dt*g*(sin(roll(i))*cos(X(5)) - cos(roll(i))*sin(pitch(i))*sin(X(5)));...
           0, 0,  0,  1, dt*g*(sin(roll(i))*sin(X(5)) + cos(roll(i))*sin(pitch(i))*cos(X(5)));...
           0, 0,  0,  0, 1];
       
    gamma = [0,0,0;...
             0,0,0;...
             sX*sqrt(dt),0,0;...
             0,sY*sqrt(dt),0;...
             0,0,sYaw*sqrt(dt)];
             
    P_Pred = Phi*P*Phi' + gamma*Q*gamma';
    
    K = P_Pred*H'*inv(H*P_Pred*H'+R);
        
    X_upd = X_Pred + K*(Y(:,i)-H*X_Pred);
    P_upd = (eye(5,5) - K*H)*P_Pred;
    
    Xm = [Xm X_upd(1,1)];
    Ym = [Ym X_upd(2,1)];
    vx = [vx X_upd(3,1)];
    vy = [vy X_upd(4,1)];
    Yaw_est = [Yaw_est X_upd(5,1)];
    
    
    X = X_upd;
    P = P_upd;    
 end

 t = 1:length(xm);

figure(1);
subplot(2,1,1);
plot(t,xm,'r-',t,Xm,'g-');
% plot(t,xm,'r-','LineWidth',1.5);
ylabel('X Co-ordinate in meters');
xlabel('Discrete time intervals k in seconds');
title('Plot of Co-ordinate measurements and Estimates');
legend('measurement', 'Estimate','Location','northeast');

subplot(2,1,2);
plot(t,ym,'r-',t,Ym,'g-');
% plot(t,ym,'r-','LineWidth',1.5);
ylabel('Y Co-ordinate in meters');
xlabel('Discrete time intervals k in seconds');
% title('Plot of Y co-ordinate measurements and estimates with time');
legend('measurement', 'Estimate','Location','northeast');
%

figure(2);
plot(t,yaw,'r-',t,Yaw_est,'g-');
% plot(t,yaw,'r-','LineWidth',1.5);
ylabel('Magnitude of Yaw Position measurements in radians');
xlabel('Discrete time intervals k in seconds');
title('Plot of Yaw measurements and estimates with time');
legend('Yaw measurement', 'Estimate','Location','northeast');