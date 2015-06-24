clear all; close all

% if a folder 'PLOTS' doesn't exist make one
if(~isequal(exist('Plots', 'dir'),7))
mkdir('Plots')
end

% load data created by run_methanogenesis
load saved_data

% PLOT OPTIONS
leg_text    = {'k', 'nup', 'chi', 'Y', 'Kac', 'm'};
x0leg_text  = {'x0(1)', 'x0(2)', 'x0(3)', 'x0(4)'};
fonsiz      = 20;
linwid      = 2;
opt_plot    = 0;


%% Morris (standard)
figure
    plot(tspan,mnt0(5:10,:),'LineWidth',linwid)
    title('Standard Morris Method','FontSize',fonsiz)
    xlabel('Time (s)','FontSize',fonsiz)
    ylabel('Mean','FontSize',fonsiz)
    set(gca,'FontSize',fonsiz)
    legend(leg_text,'FontSize',fonsiz)
    
if(opt_plot)
print -depsc Plots/figure1
end
%% Morris (no normalisation)
figure
    plot(tspan,mnt1(5:10,:),'LineWidth',linwid)
    title('Morris No Normalisation','FontSize',fonsiz)
    xlabel('Time (s)','FontSize',fonsiz)
    ylabel('Mean','FontSize',fonsiz)
    set(gca,'FontSize',fonsiz)
    legend(leg_text,'FontSize',fonsiz)

if(opt_plot)
print -depsc Plots/figure2
end
%% Morris (RSF normalization) 
figure
    plot(tspan,mnt2(5:10,:),'LineWidth',linwid)
    title('Morris RSF Normalisation','FontSize',fonsiz)
    xlabel('Time (s)','FontSize',fonsiz)
    ylabel('Mean','FontSize',fonsiz)
    set(gca,'FontSize',fonsiz)
    legend(leg_text,'FontSize',fonsiz)
    
if(opt_plot)
print -depsc Plots/figure3
end

%% TSF

figure
    plot(t,y(:,5:10),'LineWidth',linwid)
    title('TSF','FontSize',20)
    xlabel('Time (s)','FontSize',20)
    ylabel('Mean','FontSize',fonsiz)
    set(gca,'FontSize',fonsiz)
    legend(leg_text,'FontSize',20)
    
if(opt_plot)
print -depsc Plots/figure4
end
%% Comparison of the two methods without normalisation.

figure,
    plot(tspan,mnt1(5:10,:),t,y(:,5:10),'black.','LineWidth',linwid)
    title('Morris vs TSF','FontSize',20)
    xlabel('Time (s)','FontSize',20)
    ylabel('Mean','FontSize',fonsiz)
    set(gca,'FontSize',14)
    legend(leg_text,'FontSize',fonsiz)
    
if(opt_plot)
print -depsc Plots/figure5
end


%% Comparison for x0: Initial conditions

figure
    plot(tspan,mnt1(11:14,:)',t,y(:,29:32),'black.','LineWidth',linwid)
    title('x0 Morris vs TSF (before fixing)','FontSize',fonsiz)
    xlabel('Time (s)','FontSize',fonsiz)
    ylabel('Mean','FontSize',fonsiz)
    set(gca,'FontSize',14)
    ylim([-0.2,1])
    legend(x0leg_text,'FontSize',fonsiz)
    
if(opt_plot)
print -depsc Plots/figure6
end

figure
    plot(tspan,mnt3(11:14,:)',t,y(:,29:32),'black.','LineWidth',linwid)
    title('x0 Morris vs TSF (after fixing)','FontSize',fonsiz)
    xlabel('Time (s)','FontSize',fonsiz)
    ylabel('Mean','FontSize',fonsiz)
    set(gca,'FontSize',14)
    ylim([-0.2,1])
    legend(x0leg_text,'FontSize',fonsiz)

if(opt_plot)
print -depsc Plots/figure7
end