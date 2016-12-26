clear all
close all
clc

dt = 0.1;
R_full = [];
R_E = [0];
figure;
% Generate a 100 trajectories of the random walk starting with zero initial conditions 
for i = 1:100
   X = [0];
   Y = [0]; 
   R = [0];
  % Generate a random walk for 100 time steps from 0 to 10seconds
for t = 0:dt:10
    X_update = X(end)+normrnd(0,sqrt(dt));
    Y_update = Y(end)+normrnd(0,sqrt(dt));
    R_update = sqrt(X_update*X_update + Y_update*Y_update);
    X = [X X_update];
    Y = [Y Y_update];    
    R = [R R_update];
end
R_full = [R_full;R];
%Plot each trajectory of the random walk on the same figure
hold on;
plot(X,Y);
ylabel('Magnitude of X component of the Random walks');
xlabel('Magnitude of Y component of the Random walks');
title('Plot of 2D Random walk for a 100 trajectories');
end


%Problem 3 b)
% Calculate the mean of R{KdT} for each time step
R_m = [0];
for i = 1:101
    R_m = [R_m mean(R_full(:,i))];
end
% Plot all the 100 random walks along with the mean of them for each time
% step
t = 0:101;
figure;
for i = 1:100
    hold on
    plot(t,R_m,'b-',t,R_full(i,:),'g-','linewidth',1);
    ylabel('Magnitude of r(K\DeltaT)');
    xlabel('Discrete time intervals t');
    title('Plot of  r(K\DeltaT) for 100 trajectories and the mean of r(K\DeltaT)');
    legend('mean of r(K\DeltaT)', '100 trajectories of r(K\DeltaT)','Location','northeast');
end
