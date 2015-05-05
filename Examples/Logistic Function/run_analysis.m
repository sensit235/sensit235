%% run_analysis
%
% processes the files saved by run_morris and presents results.

% run_morris

% construct sensitivity measures

% find files holding experiments
file_list_es=dir('experiment-*');
nf=length(file_list_es);

% find files holding results
file_list_rs=dir('results-*');

% create array of experiments and results to feed to Process_Results

es = cell(1,nf);
rs = cell(1,nf);
% est = cell(size(r,2),nf);
% rst = cell(size(r,2),nf);
est={};
rst={};


% loop through experiments
for i=1:nf

    eval(['load ' file_list_es(i).name]);
    eval(['load ' file_list_rs(i).name]);

    % loop through times
    for j=1:size(r,2)
        est(j,i) = {A};
        rst(j,i) = {r(:,j)'};

    end

%     est(j,:) = es;
%     rst(j,:) = rs;

end

% plot a few time points
es=[est(2,:)];
rs=[rst(2,:)];
figure
plot_results
% title('short times')
% es=[est(end,:)];
% rs=[rst(end,:)];
% figure
% plot_results
% title('long times')

% tspan = linspace(0,20*24*60*60,100);
mnt = zeros(length(mn),size(est,2));
sdt = mnt;
for i=1:size(est,1)
    es=[est(i,:)];
    rs=[rst(i,:)];
    [mnt(:,i) sdt(:,i)] = Process_Results(es,rs);
end
figure
subplot(3,1,1)
plot(tspan/60/60/24,abs(mnt'))
legend({'$a$','$b$'},'Interpreter','LaTex','FontSize',20)
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
axis([0 tspan(end)/60/60/24 0 max(max([mnt sdt]))*1.01])
set(gca,'FontSize',14)

subplot(3,1,2)
plot(tspan/60/60/24,sdt')
legend({'$a$','$b$'},'Interpreter','LaTex','FontSize',20)
title('Standard Deviation vs Time','Interpreter','LaTex','FontSize',20)
axis([0 tspan(end)/60/60/24 0 max(max([mnt sdt]))*1.01])
set(gca,'FontSize',14)

% plot against solution
subplot(3,1,3)
mvt = zeros(size(est,1),1);
svt = mvt;
for i=1:size(est,1)
    mvt(i) = mean([rst{i,:}]);
    svt(i) = std([rst{i,:}]);
end
plot(tspan/60/60/24,mvt)
title('mac vs time','Interpreter','LaTex','FontSize',20)
xlabel('time (days)','Interpreter','LaTex','FontSize',20)
set(gca,'FontSize',14)

set(gcf,'PaperUnits', 'Inches', 'PaperPosition', [0 0 10 10])
print -depsc figure1


% % time averaged plot
% % cycle through parameters
% tspan = linspace(0,20*24*60*60,1000);
% for i=1:length(mn)
%     
%     mn(i) = trapz(tspan,abs(mnt(i,:)));
%     sd(i) = trapz(tspan,abs(sdt(i,:)));
%     
% end
% 
% dx=0.01; % offset to plot labels
% 
% figure
% plot(abs(mn),sd,'s','LineWidth',2)
% xlabel('Mean')
% ylabel('Standard Deviation')
% text(abs(mn)+dx.*norm(mn),sd+dx*norm(sd),{'$k$','$\nu p$','$\chi$','$Y$','$K_{ac}$','$m$'},'Interpreter','LaTex','FontSize',20)
% hold on
% plot([0 max(mn)],[0 max(mn)],'r-')
% title('Time Averaged Sensitivity Coefficients')


% % save results
% time=tspan/60/60/24;
% param_order={'$k$','$\nu p$','$\chi$','$Y$','$K_{ac}$','$m$'};
% save results-param-v-time.mat mnt sdt mvt time param_order

%% renormalise

% compute normalisation inversion values
a = [true true];

% cycle through all experiments and times results applying normalisation

ee = est;
for j=1:size(est,2) % experiments
    for i=1:size(est,1) % times
        % remove normalisations
        for k=1:size(est{i,j},1)
            pc=scale_parameters(theta,est{i,j}(k,:),a);
            est{i,j}(k,:) = pc;
        end
        
        % apply new normalisation
        [ee{i,j} sd] = Process_Results(est(i,j),rst(i,j));
        [l m n] = find(diff(est{i,j}));
        for k=1:size(ee{i,j},1)
            ee{i,j}(k) = ee{i,j}(k).*est{i,j}(l(k),k)./rst{i,j}(l(k));
        end
    end
end

% plot without normalisation
for i=1:size(est,1)
    es=[est(i,:)];
    rs=[rst(i,:)];
    [mnt(:,i) sdt(:,i)] = Process_Results(es,rs);
end
figure
% subplot(3,1,1)
plot(tspan,mnt')
legend({'$a$','$b$'},'Interpreter','LaTex','FontSize',20)
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
% axis([0 tspan(end)/60/60/24 0 max(max([mnt sdt]))*1.01])
set(gca,'FontSize',14)

% cycle through and take mean and sd of all normalised main effects


for i=1:size(est,1) % times
    mnt(:,i) = mean(cell2mat(ee(i,:)),2);
    sdt(:,i) = std(cell2mat(ee(i,:)),0,2);
end

% plot
figure
% subplot(3,1,2)
plot(tspan,mnt')
legend({'$a$','$b$'},'Interpreter','LaTex','FontSize',20)
title('Mean vs Time','Interpreter','LaTex','FontSize',20)
% axis([0 tspan(end)/60/60/24 0 max(max([mnt sdt]))*1.01])
set(gca,'FontSize',14)

% subplot(3,1,3)
% plot(tspan,sdt')
% legend({'$a$','$b$'},'Interpreter','LaTex','FontSize',20)
% title('Standard Deviation vs Time','Interpreter','LaTex','FontSize',20)
% % axis([0 tspan(end)/60/60/24 0 max(max([mnt sdt]))*1.01])
% set(gca,'FontSize',14)