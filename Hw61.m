function main
close all
clear all
clc

t = 1:10;
Z = [9055,8560, 7963, 7467, 7000, 6378, 5885, 5400, 4928, 4503];
X0 = [10000, -500, 6e7, 50,0,0,0,200,0,0,0,2e12];
r1 =1000;
r2 =500;
X_est = [];
X_upd = [];
X_meas =[]; 
P_11 = [50];
P_22 = [200];
P_33 = [2e12];
P_12 = [0];
P_13 = [0];
P_23 = [0];
R =[5,0,0;0,0,0;0,0,0];
tRec(1) = 0;
X = [10000;-500;6e7];
xRec(:,1)=X;


%discrete time steps
for dt = 1:length(Z)
    %solve ode funtion to get the discrete values at the end of time span 
    [tSol,xsol] = ode45(@radar,[dt-1 dt],X0);
    tRec=[tRec; tSol(2:end)];
    xRec=[xRec, xsol(2:end,1:3)'];
    %prediction step
    X_pred = [xsol(end,1:3)]';
    P1 = [xsol(end,4:6)];
    P2 = [xsol(end,7:9)];
    P3 = [xsol(end,10:12)];
    P_pred = [P1; P2 ;P3];
    x = X_pred(1,1);
    dx = X_pred(2,1);
    beta = X_pred(3,1);    

    
    %update step: calculate kalman gain
    h = sqrt((r1^2)+((x-r2)^2));
    % jacobian(sqrt(r1^2+(x1-r2)^2),[x1 x2 x3])
    H = [(x -r2)/( r1^2 + (x - r2)^2 )^(1/2), 0, 0];
    K = P_pred*H'*inv(H*P_pred*H'+ 5);
    
    %update step
    X_upd = X_pred + K*(Z(:,dt)-h);
    P_upd = (eye(3,3) - K*H)*P_pred;

    %store estimate values
    X_est = [X_est X_upd];
    P_11 = [P_11 xsol(2:end,4)'];
    P_22 = [P_22 xsol(2:end,8)'];
    P_33 = [P_33 xsol(2:end,12)'];
    P_12 = [P_12 xsol(2:end,5)'];
    P_13 = [P_13 xsol(2:end,6)'];
    P_23 = [P_23 xsol(2:end,9)'];
    
    
    %reset initial conditions
    X0 = [X_upd' P_upd(1,:) P_upd(2,:) P_upd(3,:)];
    clear tsol;
    clear xsol;
end
%calculate the measurement values
for i = 1:10
    X_meas(i) = (sqrt(Z(i)^2 - (r1^2)) + r2 );
end

%plots
plot([1:10],X_meas,'ro',[1:10],X_meas,'b-',tRec,xRec(1,:),'g-');
ylabel('Magnitude of Position estimates and measurements');
xlabel('Discrete time intervals t');
title('Plot of Position measurements and estimate with time');
legend('Measurement','Update',' Estimate ','Location','southwest');


figure(2);
plot(tRec,xRec(2,:));
title('Plot of Velocity estimate with time');
ylabel('Magnitude of Velocity estimate');
xlabel('Discrete time intervals t');

figure(3);
plot(tRec,xRec(3,:));
title('Plot of ballistic coefficient estimate with time');
ylabel('Magnitude of ballistic coefficient estimate');
xlabel('Discrete time intervals t');

figure(4);
plot(tRec,P_11);
title('Plot of Position variance with time');
ylabel('Magnitude of Position variance');
xlabel('Discrete time intervals t');

figure(5);
plot(tRec,P_22);
title('Plot of Velocity variance with time');
ylabel('Magnitude of  Velocity variance');
xlabel('Discrete time intervals t');

figure(6);
plot(tRec,P_33);
title('Plot of ballistic coefficient variance with time');
ylabel('Magnitude of ballistic coefficient variance');
xlabel('Discrete time intervals t');

figure(7);
plot(tRec,P_13);
title('Plot of Position-ballistic coefficient covariance with time');
ylabel('Magnitude of Position-ballistic coefficient covariance');
xlabel('Discrete time intervals t');

figure(8);
plot(tRec,P_23);
title('Plot of Velocity-ballistic coefficient covariance with time');
ylabel('Magnitude of Velocity-ballistic coefficient covariance');
xlabel('Discrete time intervals t');

figure(9);
plot(tRec,P_12);
title('Plot of Position-Velocity covariance with time');
ylabel('Magnitude of Position-Velocity covariance');
xlabel('Discrete time intervals t');

end

%ode function
function dx = radar(t,X)
b = [0;0;sqrt(1000)];
X1 = X(1,1);X2 = X(2,1);X3 = X(3,1);
P11 = X(4,1);P12 = X(5,1);P13 = X(6,1);
P21 = X(7,1);P22 = X(8,1);P23 = X(9,1);
P31 = X(10,1);P32 = X(11,1);P33 = X(12,1);

dotX1 =X2;
dotX2 = ((1220*exp(-X1/10263)*X2*X2)/(2*X3))- 9.8;
dotX3 = 0;

P = [P11,P12,P13;P21,P22,P23;P31,P32,P33];
% jacobian([x2;((1220*exp(-x1/10263)*x2^2)/(2*x3))- 9.8;0],[x1 x2 x3])
A = [0,1,0;-(610*X2^2*exp(-X1/10263))/(10263*X3), (1220*X2*exp(-X1/10263))/X3, -(610*X2^2*exp(-X1/10263))/X3^2;0,0,0];
dotP  = A*P + P*A' + b*b';

dx(1,1) = dotX1;dx(2,1) = dotX2;dx(3,1) = 0;
dx(4,1) = dotP(1,1);dx(5,1) = dotP(1,2);dx(6,1) = dotP(1,3);
dx(7,1) = dotP(2,1);dx(8,1) = dotP(2,2);dx(9,1) = dotP(2,3);
dx(10,1) = dotP(3,1);dx(11,1) = dotP(3,2);dx(12,1) = dotP(3,3);
end

