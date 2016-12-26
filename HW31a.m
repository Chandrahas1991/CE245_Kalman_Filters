clear all
close all
clc

%input the parameters of the linear dynamical system
A = [1.5, 1;-0.7, 0];
B = [1;0.5];
X = [];
X(:,1) = [0;0];
%initialize the state covariance matrix
P= [0,0;0,0];
G = [1;0.5];
P_11 = 0;
P_22 = 0;
P_12 = 0;
% propogate the covariance matrix of system by 30 iterations
for k =1:30
    P = A*P*A' + G*G';
    P_11 = [P_11 P(1,1)];
    P_22 = [P_22 P(2,2)];
    P_12 = [P_12 P(1,2)];
end
%plot the state covariance transition from the initial state 
k = 1:31;
figure;
hold on
plot(k,P_11,'r-',k,P_22,'g-',k,P_12,'b-');
ylabel('Magnitude of elements Covarience matrix P');
xlabel('Discrete time intervals k');
title('Plot of State Covariance Transition');
legend('P_{11}', 'P_{22}','P_{12} = P_{21}','Location','northeast');

