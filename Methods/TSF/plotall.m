clear all;close all;clc

load res_num.mat

figure(1)
subplot(211)
plot(t_num,y_num(:,13))
title('K_{Ac}')
subplot(212)
plot(t,y(:,5))

figure(2)
subplot(211)
plot(t_num,y_num(:,14))
title('m')
subplot(212)
plot(t,y(:,6))

figure(3)
subplot(211)
plot(t_num,y_num(:,11))
title('Y')
hold on
subplot(212)
plot(t,y(:,7))

figure(4)
subplot(211)
plot(t_num,y_num(:,8))
title('nup')
hold on
subplot(212)
plot(t,y(:,8))

figure(5)
subplot(211)
plot(t_num,y_num(:,7))
hold on
title('k')
subplot(212)
plot(t,y(:,9))

figure(6)
subplot(211)
plot(t_num,y_num(:,10))
hold on
title('chi')
subplot(212)
plot(t,y(:,10))

figure(101)
subplot(211)
plot(t_num,y_num(:,1))
hold on
title('state 1')
subplot(212)
plot(t,y(:,1))

figure(102)
subplot(211)
plot(t_num,y_num(:,2))
hold on
title('state 2')
subplot(212)
plot(t,y(:,2))

figure(200)
h1=plot(t_ic,y_ic(:,1),'LineWidth',2)
hold on
h2=plot(t_ic,y_ic(:,2),'LineWidth',2)
h3=plot(t_ic,y_ic(:,3),'LineWidth',2)
h4=plot(t_ic,y_ic(:,4),'LineWidth',2)
legend([h1,h2,h3,h4], 'dm_{Ac}/dx_1','dm_{Ac}/dx_2','dm_{Ac}/dx_3','dm_{Ac}/dx_4')
title('SENSITIVITIES wrt IC')