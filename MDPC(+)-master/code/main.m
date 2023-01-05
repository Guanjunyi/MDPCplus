% Please kindly cite the paper Junyi Guan, Sheng li, Xiongxiong He, and Jiajia Chen 
%"Clustering by fast detection of main density peaks within a peak digraph" 
% Information Sciences,2023

% The code was written by Junyi Guan in 2022.
clear;close all;clc;
%% load dataset
load data/Agg
data_with_lable = Agg;
%% deduplicate data
data_x = unique(data_with_lable,'rows');
if size(data_x,1) ~= size(data_with_lable,1)
    data_with_lable = data_x;
end
lable = data_with_lable(:,end);
data = data_with_lable(:,1:end-1);
%% data preprocessing
data=(data-min(data))./(max(data)-min(data));%% data normalization
data(isnan(data))=0;
%% MDPC+ clustering
[CL,NC,centers,runtime] = MDPC_Plus(data);
%% evaluation
[AMI,ARI,FMI] = Evaluation(CL,lable);
%% show result
resultshow(data,CL,centers);
%% clustering result
result = struct;
result.NC = NC;
result.AMI = AMI;
result.ARI = ARI;
result.FMI = FMI;
result.runtime = runtime;
result