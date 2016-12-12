function main
clear all;
close all
clc;
% Set the Upper Bound of the time interval
T = 500;
t = [];
t(1) = 0;
k = 1;
% Generate a random exponentialy varying time interval for the disturbance 
    while(t(k)< T)
        k = k + 1;
        % The parameter for exponential distribution is the mean
        t(k) = t(k-1) + random('exp',50);
    end
% Set the last 
t(k) = T;
Num_Events = k;
% Set inital conditions of the system
x_o = [0 0];
% Initialize empty arrays to store the solutions
xi = [];
ti =[];
% Calculate the solutions for each time interval 
    for k = 1:(Num_Events - 1)
        [ts,xs] = ode45(@deq,[t(k): 0.5 : t(k+1)],x_o);
        ti = cat(1,ti,ts);
        xi = cat(1,xi,xs);
        ts = [];
        xs = [];
        % Add the normally distributed random event('Kick') to the initial
        % condition
        f =random('normal',0,0.2); 
        x_o = xi(end) + [0,f];
end

% Plot the solutions against time.
plot(ti,xi(:,1),'k-','Linewidth',1)
xlabel('Time, t (sec)'),ylabel('Angular Position, x (rad)'),grid on;
title('Plot of Angular Position Vs Time');
figure;
plot(ti,xi(:,2),'r-','Linewidth',1);
xlabel('Time, t (sec)'),ylabel('Angular Velocity, x (rad)'),grid on;
title('Plot of Angular Velocity Vs Time');
figure;

%Plot the Angular Position against Angular Velocity.
plot(xi(:,1),xi(:,2),'g-','linewidth',1);
xlabel('Angular Position, x (rad)'),ylabel('Angular Velocity, x (rad)'),grid on;
title('Plot of Angular Position Vs Angular Velocity');

% Define a function for the system of Differential Equations
    function dxdt = deq(t,x)
        dxdt = [x(2); -0.1*sin(x(1))];
    end

end


 