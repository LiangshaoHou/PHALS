clc; close all;clear all;
addpath(genpath(pwd));

% Generate similarity matrix: 
% load('ORL_64x64.mat'); % load('Coil20.mat');% Coil20.mat data set
% A=Image_Similarity(fea);

load('ORL_Similar.mat');% load similarity matrix and label of ORL dataset
% load('COIL_Similar.mat'); % load similarity matrix and label of COIL dataset
n=size(A,1);
label=gnd;
r=max(label);


% Initialization
U0 = rand(n,r);
H = sqrt(max(trace(U0'*A*U0),0))/norm(U0'*U0, 'fro')*U0;
U0 = abs(H);

% Test PHALS
Maxtime = 360; % max running time
Maxiter=1000;
tic
[U,Obj_PHALS,gap_PHALS,time_PHALS]=PHALS(A,U0,Maxiter, Maxtime);
toc

figure(1)
plot(time_PHALS,Obj_PHALS);
xlabel('CPU');
ylabel('Relative Error');

figure(2)
semilogy(time_PHALS,gap_PHALS);
xlabel('CPU');
ylabel('Optimality Gap');

% If ground truth labels are provided, update accuracy
[~,plabel]=max(U,[],2); % Compute predicted labels from U
temp=ClusteringMeasure(label,plabel); % Measure clustering accuracy
Acc=temp(1); % Store accuracy


