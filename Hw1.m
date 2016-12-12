function main
% Run Two Cases with Different Initial Positions
clc;clear all;clf
g=9.81;m=100;L=1;Io=50;
K=m*g*L/Io;
xo=[pi/3,0];
[t,x]=ode45(@DE2,[0:0.01:5],xo);
plot(t,x(:,1),'k--','Linewidth',1.5)
xlabel('Time(t)'),ylabel('Angular Motion, \theta (rad)'),grid on
title('Numerical Solution of d^2\theta/dt^2+Ksin\theta=0')
axis([0,5,-1.5,2])
hold on
xo=[pi/10,0];
[t,x]=ode45(@DE2,[0:0.01:5],xo);
plot(t,x(:,1),'k-','Linewidth',1.5)
% Plot Linearized Solution of Previous Two Cases
to=pi/3,td=0;
T=0:0.05:5;
th=td*sin(sqrt(K)*T)+to*cos(sqrt(K)*T);
plot(T,th,'r*','Linewidth',1)
to=pi/10,td=0;
T=0:0.2:5;
th=td*sin(sqrt(K)*T)+to*cos(sqrt(K)*T);
plot(T,th,'b*','Linewidth',1)
legend('\theta(0) = \pi/3','\theta(0) = \pi/10'...
,'Linear \theta(0) = \pi/3','Linear \theta(0) = \pi/10')
function dxdt=DE2(t,x)
g=9.81;m=100;L=1;Io=50;
K=m*g*L/Io;
dxdt=[x(2);-K*sin(x(1))];