

% We assume that there are N=3 network nodes X1,X2, and
% We assume that there are N=3 network nodes X1,X2, and
% We assume that there are N=3 network nodes X1,X2, and X3
% in a network with two edges:
% X1 -> X2
% X1 -> X3
% We have sampled (X1,X2,X3) for 500 time points.
% There are two transition times:
% t1 = 200 for X2, and t2 =300 for X3.

% Install HMDBN
addpath('.\HMDBN') 
% Load testing data
data=load('data.txt'); 
% Run HMDBN and calculate the running time
time1 = cputime;
[dag hiddenGraph_Ps SampleDistribution] = hmdbn_structEM(data);
ElapsedTime = cputime - time1;

% Demonstrate the distribution of each edge over time.
hmdbn_demonstrate(data,hiddenGraph_Ps,SampleDistribution);
