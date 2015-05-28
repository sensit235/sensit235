function arun_sens_logistic_num

addpath ../../External/DERIVESTsuite/

% Simulation interval
tspan = [0,18];
% Parameters 
a  = 0.8; b  = 0.1;
pars = [a b];
% Initial condition
x0 = 0.3;       
xinit = [x0 0 0];
% Solution
[t0,y0] = ode15s(@rhs_sens_logistic_num,tspan,xinit,[],pars);

% Analytic
[t1,y1] = ode15s(@rhs_sens_logistic,tspan,xinit,[],pars);


% Plot
figure(1)
h1=plot(t0,y0(:,2),'LineWidth',2);
hold on
h2=plot(t0,y0(:,3),'LineWidth',2);
grid on
set(gca,'FontSize',14)
title('Sensitivities')
h3=plot(t1,y1(:,2),'x','LineWidth',2);
h4=plot(t1,y1(:,3),'x','LineWidth',2);
legend([h1 h2 h3 h4],'dx/da','dx/db','analytic','analytic')

end


function dy = rhs_sens_logistic(~,v,pars)
% pars
a = pars(1);
b = pars(2);
% define sensitivities
s11=v(2);
s12=v(3);
dy = zeros(3,1);

% Model: xd = F(x)
dy(1) = a*v(1) - b*v(1)^2;
% Sensitivity equation 
dy(2) =  v(1) + s11*(a - 2*b*v(1));
dy(3) = -v(1)^2 + s12*(a - 2*b*v(1));
end

% numeric version
function dy = rhs_sens_logistic_num(~,v,pars)

% define sensitivities
dy = zeros(3,1);

% F(x)
s = @(x,p) p(1).*x - p(2).*x.^2;

% Model
dy(1) = s(v(1),pars);

% dF/dx
dfdx = derivest(@(x) s(x,pars),v(1));

% dF/da
dfda = jacobianest(@(x) s(v(1),x),pars);

% Sensitivity equation
dy(2) = dfdx*v(2) + dfda(1)
dy(3) = dfdx*v(3) + dfda(2)

end
