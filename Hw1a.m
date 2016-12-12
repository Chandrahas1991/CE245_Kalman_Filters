function main 
clear all;
clc;
clf;
%set initial conditions.
xo = [pi/6,0];

% Solve the system of equations using the ode45 function.
[t,x] = ode45(@deq,[0 500],xo);

% Plot the solutions against time.
plot(t,x(:,1),'k-',t,x(:,2),'r-','Linewidth',1);
legend('Angular Position', 'Angular Velocity','Location','northeast');
xlabel('Time (t)'),ylabel('Angular Motion, x (rad)'),grid on;
title('Numerical Solution of d^{2}x/dt^{2}+Ksin(x)=0');
figure;

%Plot the Angular Position against Angular Velocity.
plot(x(:,1),x(:,2),'k-','linewidth',1);
xlabel('Angular Position X1 (rad)'),ylabel('Angular Velocity X2 (rad/sec)'),grid on;
title('Plot of Angular Position X1 (rad) VS Angular Velocity X2 (rad/sec)');
end

% Define a function for the system of Differential Equations
function dxdt = deq(t,x)
dxdt = [x(2); -0.1*sin(x(1))];
end

