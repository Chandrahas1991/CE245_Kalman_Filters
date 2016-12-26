clear all
close all
clc

%input the parameters of the linear dynamical system
A = [1.5, 1;-0.7, 0];
B = [1;0.5];
X = [];
X(:,1) = [0;0];
%problem 1 b)
% propogate the system by for 5 trajectories each for 500 discrete time steps 
figure;
for i = 1:5
    for K = 2:500
        X(:,K) = A*X(:,K-1) + B*random('normal',0,1);
    end
k = 1:500;
hold on
subplot(2,1,1),plot(k,X(1,:));
ylabel('Magnitude of x_{1}(k)');
xlabel('Discrete time intervals k');
hold on
subplot(2,1,2),plot(k,X(2,:));
ylabel('Magnitude of x_{2}(k)');
xlabel('Discrete time intervals k');
title('State propogation with discrete time interval k');
end

%problem 1 c)
%Generate 50 trajectories starting with zero initial conditions and each for 500 discrete time steps
Z(:,1) = zeros(100,1);   
for i = 1:2:50
    for K = 2:500
        Z(i:(i+1),K) = A*Z(i:(i+1),K-1) + B*random('normal',0,1);
    end
end
%store the values of x(1) and x(2) at the 500th time step from all the trajectories and use them to calculate the standard deviation. 
X1_merge = [];
X2_merge = [];
for i = 1:2:100
    X1_merge = [X1_merge Z(i,500)];
end
for i = 2:2:100
    X2_merge = [X2_merge Z(i,500)];
end
% Standard deviations of X_1(500) and X_2(500)
SD_X1 = sqrt(var(X1_merge));
SD_X2 = sqrt(var(X2_merge));
% SD of x(1) ~ 3.5 and SD of x(2) ~ 2.5

%problem 1 d)
P= [ 1, 0.5; 0.5, 0.25];
G = [1;0.5];
for i =1:50
    P = A*P*A' + G*G';
end
SD_X1_new = sqrt(P(1,1));
SD_X2_new = sqrt(P(2,2));
