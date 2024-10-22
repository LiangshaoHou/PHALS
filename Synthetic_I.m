%% Comparing our algorithm with SymHALS and SymANLS
clc; close all;clear all;
addpath(genpath(pwd));
n=1000;
r=50;


% Synthetic dataset
A=abs(randn(n,r));
A=A*A';

% Initialization
U0 = rand(n,r);
H = sqrt(max(trace(U0'*A*U0),0))/norm(U0'*U0, 'fro')*U0;
U0 = abs(H);

% Runing PHALS
% Test PHALS
Maxtime = 360; % max running time
Maxiter=300; 
[U,Obj_Val,Opt_Gap,time]=PHALS(A,U0,Maxiter, Maxtime);