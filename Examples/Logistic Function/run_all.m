%% script for running logistic example

close all
clear all

addpath ../../Methods/Morris' Method'/

run_gen_exp(4,4,4);

a = [true true true true];
leg_text = {'$a$','$b$','$c$','$d$'};

tspan = linspace(0,16,100); % time span or vector of times
theta = [0.8,0.1,0.3,1.0]; % parameters
run_morris(tspan,theta,a)

run_analysis