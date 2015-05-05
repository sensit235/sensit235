%% plot_results
%
% If this script is called after Process Results the mean and standard
% deviation can be used to plot the results.

[mn sd] = Process_Results(es,rs);

dx=0.01; % offset to plot labels

plot(abs(mn),sd,'s','LineWidth',2)
xlabel('Mean')
ylabel('Standard Deviation')
% text(abs(mn)+dx.*norm(mn),sd+dx*norm(sd),cellstr(num2str([1:length(mn)]')))
% text(abs(mn)+dx.*norm(mn),sd+dx*norm(sd),{'$k$','$\nu p$','$\chi$','$Y$','$K_{ac}$','$m$'},'Interpreter','LaTex','FontSize',20)
