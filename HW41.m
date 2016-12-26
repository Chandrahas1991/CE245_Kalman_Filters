close all
clear all
clc
% Initializations of state space 
X1(1) = 10;
X2(1) = 20;
XN = [10;20];
% Initializations of expected values
EXN = [10;20];
U = [0;1];
XN_full = XN;
EXN_full = XN;
% Initializations of variances 
P = diag([40,10000]);
PX1 = [P(1,1)];
PX2 = [P(2,2)];
for K = 1:14
    W1 = normrnd(0,1);
    W2 = normrnd(0,1);
    W = [W1;W2];
    Q = W*W';
    O = [0.2, 0.4; -0.4,1];
    % State space equations
    XN = O*XN + U + W;
    XN_full = [XN_full XN];
    %Expected values equations
    EXN = O*EXN + [0;1];
    EXN_full = [EXN_full EXN];
    % Variance equations
    P = O*P*(O') + Q;
    PX1 = [PX1 P(1,1)];
    PX2 = [PX2 P(2,2)];
    
end

K = 0:14;

figure;
plot(K,EXN_full(1,:),K,EXN_full(2,:));
legend('X_{1}(K)', 'X_{2}(K)','Location','northeast');
ylabel('Magnitude of Expected values of state variables X_{1} and X_{2}');
xlabel('Discrete time intervals K');
title('Plot of state variables transistion');

figure;
plot(K,PX1,K,PX2);
legend('P_{11}', 'P_{22}','Location','northeast');
ylabel('Magnitude of covariances of state variables X_{1} and X_{2}');
xlabel('Discrete time intervals K');
title('Transistion of covariance of state variables');
