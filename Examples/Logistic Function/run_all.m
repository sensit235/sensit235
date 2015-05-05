%% script for running logistic example

close all
clear all

addpath ../../Methods/Morris' Method'/

run_gen_exp(2,10,2);

tspan = linspace(0,16,100); % time span or vector of times
theta = [0.8,0.1]; % parameters
initial = [0.1]; % initial condition
run_morris(tspan,theta,initial)

run_analysis