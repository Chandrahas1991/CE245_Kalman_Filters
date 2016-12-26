function main 
close all
clear all
clc
%Initialize the values
S_0=100;
V_0=100;
[t CovP]=ode45(@CovMat,[0 3],[S_0;0;V_0]);
P_11=CovP(:,1);
P_12=CovP(:,2);
P_22=CovP(:,3);
plot(t,P_11);
hold on
plot(t,2*S_0*ones(length(t)),'r-');
xlabel('Time t');
ylabel('Transistion of state variance var \{s(t)\} = P_{11}');
figure(2)
plot(t,P_22)
xlabel('Time t');
ylabel('Transistion of state variance var \{v(t)\} = P_{22}');
figure(3)
plot(t,P_12)
xlabel('Time t');
ylabel('Transistion of state covariance cov \{s(t)v(t)\} = P_{12}= P_{21}');
end
% ODE function solver.
function dCovX=CovMat(t,X)
  P(1,1)=X(1,1);
  P(1,2)=X(2,1);
  P(2,1)=P(1,2);
  P(2,2)=X(3,1);
  A=[0 1; 0 -2];
  W=[0; 20];
  dP=A*P+P*A'+W*W';
  dCovX=[dP(1,1);dP(1,2);dP(2,2)];
end
