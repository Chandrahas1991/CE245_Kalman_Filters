close all
clear all
clc

res_hist =[]; 
D3D = load('Data3D.mat');
D1 = load('dataset1.mat');
D2 = load('dataset2.mat');
D3 = load('dataset3.mat');
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
thrust = D3D.thrust;
s = 255/6000;
g =9.81;
sX = 0.1;
sY = 0.1;
sZ = 0.1;
sXd = 2*sX/(dtrec(1)^2);
sYd = 2*sY/(dtrec(1)^2);
sZd = 2*sZ/(dtrec(1)^2);
sYaw = degtorad(0.01);

for c0 = -199:2:199
   
c = (0.000409*(s^2) + 0.1405*s -0.099)/c0;
sum_res = 0;

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
Yaw_est = [yaw(1)];

P = [sX,0,0,0,0,0,0;
     0,sY,0,0,0,0,0;
     0,0,sZ,0,0,0,0;
     0,0,0,sXd,0,0,0;
     0,0,0,0,sYd,0,0;
     0,0,0,0,0,sZd,0;
     0,0,0,0',0,0,sYaw];
 
for i = 2:length(dtrec)
    X(1) = X(1)+ dtrec(i)*X(4);
    X(2) = X(2)+ dtrec(i)*X(5);
    X(3) = X(3)+ dtrec(i)*X(6);
    X(4) = X(4)+ dtrec(i)*c*g*(sin(X(7))*sin(roll(i))+cos(X(7))*cos(roll(i))*sin(pitch(i)));
    X(5) = X(5)+ dtrec(i)*c*g*(sin(X(7))*cos(roll(i))*sin(pitch(i))-cos(X(7))*sin(roll(i)));
    X(6) = X(6)+ dtrec(i)*c*(g*cos(roll(i))*cos(pitch(i))-g);
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
      
    %calculating residuals
    
    sum_res = sum_res + ((xm(i)- X_Pred(1,1))^2 + (ym(i) - X_Pred(2,1))^2 + (zm(i) - X_Pred(3,1))^2);
       
    X = X_upd;
    P = P_upd;    
end

residual = sum_res/length(xm);
res_hist = [res_hist residual];
end

%Plot
c0 = -199:2:199
plot(c0,res_hist);
xlabel('C0');
ylabel('Sum of square of residuals in meters^{2}');

