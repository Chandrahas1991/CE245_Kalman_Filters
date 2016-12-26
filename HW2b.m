clear all
close all
clc
%create two uniform probablity ditribution functions with a lower bound of
%-0.25 and an upper bound of 0.25
P_n = makedist('uniform','lower',-0.25,'upper',0.25);
P_h = makedist('uniform','lower',-0.25,'upper',0.25);
% Define the number of sequences.
k = 100000;
%sample both the probability distributions 'K'times. 
samples1 = random(P_n,k,1);
samples2 = random(P_h,k,1);
%Add both of these samples to obtain a set of samples for the random
%variable mk.
conv = samples1 + samples2;
% Display the data obtained as an histogram which show that the convoultion of
%the two uniform PDF results in a triangular PDF.
[f x] = hist(conv,100);
% Add two corner conditions(p(mk) = 0 when mk>0.5 and mk< -0.5) to close 
% the probability distribution fuction which now resembles the probablity
% distribution of the given random variable mk
f = [0 f 0];
x = [-0.51 x 0.51];
figure;
% Normalize the graph such that area under the graph is 1,
% and plot the PDF as both histogram and line graph
bar(x,f/trapz(x,f));
ylabel('Probability Density');
xlabel('Random Variable m_{k}')
title('Histogram of the PDF corresponding to the random variable m_{k}')
figure;
plot(x,f/trapz(x,f));
ylabel('Probability Density');
xlabel('Random Variable m_{k}')
title('Probability distribution function of random variable m_{k}')
% The mean of the PDF is given by, (mean_mk ~ 0)
mean_mk = mean(conv);
% The varience of the PDF is given by, (var_mk ~ 0.04 = 1/24)
var_mk = var(conv);
% The varience of state of the stochastic system after 'k' sequences is
% given by, (for k =100000, var_xk ~ 4170)
var_xk = k*var_mk;