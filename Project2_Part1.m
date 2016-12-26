close all
clear all
clc

D3D = load('Data3D.mat');
dtrec = D3D.dtrec(45:end);
pitch = degtorad(D3D.pitch(45:end));
roll = degtorad(D3D.roll(45:end));
yaw = degtorad(D3D.yawm(45:end));
xm = D3D.xm(45:end);
ym = D3D.ym(45:end);
zm = D3D.zm(45:end);
vx(1) = (xm(2)-xm(1))/dtrec(1);
vy(1) = (ym(2)-ym(1))/dtrec(1);
vz(1) = (zm(2)-zm(1))/dtrec(1);
Y = [xm';ym';zm';yaw'];
thrust = D3D.thrust(45:end);
g =9.81;
sX = 0.1;
sY = 0.1;
sZ = 0.1;
sXd = 2*sX/(dtrec(1)^2);
sYd = 2*sY/(dtrec(1)^2);
sZd = 2*sZ/(dtrec(1)^2);
sYaw = degtorad(0.01);
P =[];

Q = [1,0,0,0;...
     0,1,0,0;...
     0,0,1,0;...
     0,0,0,1];
R = [0.01^2,0,0,0;...
     0,0.01^2,0,0;...
     0,0,0.01^2,0;...
     0,0,0,(pi/180)^2];
H = [1,0,0,0,0,0,0;...
     0,1,0,0,0,0,0;...
     0,0,1,0,0,0,0;...
     0,0,0,0,0,0,1];
X  = [xm(1);ym(1);zm(1);vx(1);vy(1);vz(1);yaw(1)];
Xm = [xm(1)];
Ym = [ym(1)];
Zm = [zm(1)];
xvel_hist = [vx(1)];
yvel_hist = [vy(1)];
zvel_hist = [vz(1)];

Yaw_est = [yaw(1)];

P = [sX,0,0,0,0,0,0;
     0,sY,0,0,0,0,0;
     0,0,sZ,0,0,0,0;
     0,0,0,sXd,0,0,0;
     0,0,0,0,sYd,0,0;
     0,0,0,0,0,sZd,0;
     0,0,0,0,0,0,sYaw];
 
 
for i = 2:length(dtrec)
    X(1) = X(1)+ dtrec(i)*X(4);
    X(2) = X(2)+ dtrec(i)*X(5);
    X(3) = X(3)+ dtrec(i)*X(6);
    X(4) = X(4)+ dtrec(i)*g*(sin(X(7))*sin(roll(i))+cos(X(7))*cos(roll(i))*sin(pitch(i)));
    X(5) = X(5)+ dtrec(i)*g*(sin(X(7))*cos(roll(i))*sin(pitch(i))-cos(X(7))*sin(roll(i)));
    X(6) = X(6)+ dtrec(i)*(g*cos(roll(i))*cos(pitch(i))-g);
    X(7) = wrapToPi(X(7));
    
    X_Pred = X;
%   jacobian([X1+ dtrec(i)*X4,X2+ dtrec(i)*X5,X3+ dtrec(i)*X6,X4+ dtrec(i)*g*(sin(X7)*sin(roll(i))+cos(X7)*cos(roll(i))*sin(pitch(i))),X5+ dtrec(i)*g*(sin(X7)*cos(roll(i))*sin(pitch(i))-cos(X7)*sin(roll(i))),X6+ dtrec(i)*g*cos(roll(i))*cos(pitch(i))-g, X7],[X1 X2 X3 X4 X5 X6 X7])
 
    Phi =[ 1, 0, 0,  dtrec(i),  0,  0,  0;...
           0, 1, 0,  0,  dtrec(i),  0,  0;...
           0, 0, 1,  0,  0, dtrec(i),   0;...
           0, 0, 0,  1,  0,  0, g*dtrec(i)*(sin(roll(i))*cos(X(7)) - cos(roll(i))*sin(pitch(i))*sin(X(7)));...
           0, 0, 0,  0,  1,  0, g*dtrec(i)*(cos(roll(i))*sin(pitch(i))*cos(X(7)) + sin(roll(i))*sin(X(7)));...
           0, 0, 0,  0,  0,  1, 0;...
           0, 0, 0,  0,  0,  0, 1];
    gamma = [0,0,0,0;...
             0,0,0,0;...
             0,0,0,0;...
             sX*sqrt(dtrec(i)),0,0,0;...
             0,sY*sqrt(dtrec(i)),0,0;...
             0,0,sZ*sqrt(dtrec(i)),0;...
             0,0,0,sYaw*sqrt(dtrec(i))];
             
    P_Pred = Phi*P*Phi' + gamma*Q*gamma';
    
    K = P_Pred*H'*inv(H*P_Pred*H'+R);
        
    X_upd = X_Pred + K*(Y(:,i)-H*X_Pred);
    P_upd = (eye(7,7) - K*H)*P_Pred;
    
    Xm = [Xm X_upd(1,1)];
    Ym = [Ym X_upd(2,1)];
    Zm = [Zm X_upd(3,1)];
    vx = [vx X_upd(4,1)];
    vy = [vy X_upd(5,1)];
    vz = [vz X_upd(6,1)];
    
    xvel = (xm(i) - xm(i-1))/dtrec(i);
    yvel = (ym(i) - ym(i-1))/dtrec(i);
    zvel = (zm(i) - zm(i-1))/dtrec(i);
    
    
    xvel_hist = [xvel_hist xvel];
    yvel_hist = [yvel_hist yvel];
    zvel_hist = [zvel_hist zvel];
    
   
    Yaw_est = [Yaw_est X_upd(7,1)];
    
    
    X = X_upd;
    P = P_upd;    
end

%Plots

t = 1:length(xm);

figure(1);
subplot(3,1,1);
plot(t,xm,'r-',t,Xm,'g-');
ylabel('X Co-ordinate in meters');
xlabel('Discrete time intervals t(sec)');
title('Plot of Co-ordinate measurements and Estimates');
a = get(gca,'Children');


subplot(3,1,2);
plot(t,ym,'r-',t,Ym,'g-');
ylabel('Y Co-ordinate in meters');
xlabel('Discrete time intervals t(sec)');
b = get(gca,'Children');

subplot(3,1,3);
plot(t,zm,'r-',t,Zm,'g-');
ylabel('X Co-ordinate in meters');
xlabel('Discrete time intervals t(sec)');
c = get(gca,'Children');
o = [a;b;c];
legend(o,'measurement', 'Estimate','Orientation','horizontal');


figure(2);
plot(t,yaw,'r-',t,Yaw_est,'g-');
ylabel('Magnitude of Yaw Position measurements in radians');
xlabel('Discrete time intervals t(sec)');
title('Plot of Yaw measurements and estimates with time');
legend('Yaw Measurement', 'Yaw Estimate','Location','northeast');

figure(3);
plot3(Xm,Ym,Zm,'g-',xm,ym,zm,'r-');
ylabel('Y-Axis Co-ordinates');
xlabel('X-Axis Co-ordinates');
zlabel('Z-Axis Co-ordinates');
title('Quadcopter Trajectory');
legend('Estimate','Measurement','Location','northeast');

figure(4);
subplot(3,1,1);
plot(t,xvel_hist,'r-',t,vx,'g-');
d = get(gca,'Children');
title('Plot of Velocity Estimates in meters/sec');
ylabel('X co-ordinate Velocity');
xlabel('Discrete time intervals t(sec)');

subplot(3,1,2);
plot(t,yvel_hist,'r-',t,vy,'g-');
ylabel('Y Co-ordinate Velocity');
xlabel('Discrete time intervals t(sec)');
e = get(gca,'Children');

subplot(3,1,3);
plot(t,zvel_hist,'r-',t,vz,'g-');
ylabel('Z Co-ordinate Velocity');
xlabel('Discrete time intervals t(sec)');
f = get(gca,'Children');
u=[d;e;f];
legend(u,'Estimate','measurement', 'Orientation','horizontal');













