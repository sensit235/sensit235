function x_methane_sens_num

close all


addpath DERIVESTsuite/

% load data
t_data=[0 2 4 5 6 10 12 14 16 18]*24*60*60;
y_data=[3.68 2.54 1.90 1.85 1.61 1.31 1.28 1.17 1.02 0.97]*1e-3;

% Simulation interval
days = 18;
td=(0:0.01:1)*days*24*60*60;


% load parameters
DG0 = -15802.1961;  % standard value of Gibbs free energy at T=310.15K (J/mol)
R   = 8.3145;       % gas constant (J/mol/K)
%k   = 2.5e-6;       % (mol/g/s) **
k   = 2.5e-5;       % (mol/g/s) **

nup = 0.5;          % (per reaction) **
DGp = 45000;        % phosphorylation energy (J/(mol ATP)) **
chi = 2;            % (per reaction) **
Y   = 2.1;          % (Ymax) (g/mol) **
T   = 310.15;       % physiological temperature (K)
Kac = 5e-3;         % from Chewy  % (Kd) half saturation constant (molal) **
m   = 2.2e-7;       % from Chewy      % (D) specific maintenance rate (1/s) **

%Kn  = 0;           % effect of nutrient (molal)

% x0: fixed 
x0_mca  = 3.5e-3;
x0_hco3 = 45.2e-3; 
x0_ch4  = 1e-6; 
x0_x    = 0.009;

% params vectors
par_fis     = [DG0 R];
par_mic     = [k nup DGp chi Y T];
par_ext     = [Kac m];
% pars: fixed / estimable
par_fix = [par_fis par_mic];
par_est = par_ext;
% x0: fixed / estimable
x0_mca  = 3.5e-3;
x0_hco3 = 45.2e-3; 
x0_ch4  = 1e-6; 
x0_x    = 0.009;
% x0_fix = [];
% x0_est = [x0_mca, x0_hco3, x0_ch4, x0_x];
x0_fix = [x0_mca, x0_hco3, x0_ch4, x0_x];
x0_est = [];
%x0 = [3.5e-3 45.2e-3 1e-6 0.009];          % [mmolal mmolal mmolal g/kg]

% tspan
tspan = [td(1),td(end)];

options = odeset('RelTol',1e-6,'AbsTol',1e-6);

% Solve
x0=[x0_fix x0_est 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
par_init = [par_fix  par_est];
[t,y]   = ode15s(@rhs_sens_methane,tspan,x0,options,par_init);

% Numerical coputations
xinit2 = [x0_fix zeros(1,40)];
pp = [DG0, R, k, nup, DGp, chi, Y, T, Kac, m];
n_s = 4;
n_p = 10;
[t_num,y_num] = ode15s(@rhs_num,tspan,xinit2,options,pp,n_s,n_p);

save res_num.mat t y t_num y_num


% sen1=norm(Kac*y(:,5)./y(:,1),2)    
% sen2=norm(m*y(:,6)./y(:,1),2)    
  
end

function dy = rhs_num(t,v,pars,n_s,n_p)
%
% v     : state and sensitivity variables
% pars  : parameter values
% n_s   : number of states
% n_p   : number of parameters

%pp       = [DG0, R, k, nup, DGp, chi, Y, T, Kac, m];
% DG0     = pars(1);
% R       = pars(2);
% k       = pars(3);
% nup     = pars(4);
% DGp     = pars(5);
% chi     = pars(6);
% Y       = pars(7);
% T       = pars(8);
% Kac     = pars(9);
% m       = pars(10);

rr = @(x,p) p(3)*x(4)*x(1)/(p(9) + x(1))*...
     (1 - exp(...
     (p(1)+p(2)*p(8)*log(x(2)*x(3)/x(1))+p(4)*p(5))/...
     p(6)/p(2)/p(8)));
 
ss = @(x,p)[-rr(x,p);... 
             rr(x,p);...
             rr(x,p);...
             p(7)*rr(x,p) - p(10)*x(4)];

%ss = @model;

% Model
dy_state = ss([v(1:n_s)],pars); 

% Jacobian df/dx
    %dfdx = derivest(@(x) s(x,pars),v(1));
dfdx = jacobianest(@(x) ss(x,pars),[v(1:n_s)]);

%size(dfdx)

% dF/da
    %dfda = jacobianest(@(x) s(v(1),x),pars);
dfda = jacobianest(@(x) ss([v(1:n_s)],x),pars);


% Sensitivity equation
    %dy(5) = dfdx*v(2) + dfda(1);
    %dy(6) = dfdx*v(3) + dfda(2);
    
v40     = v(n_s+1:n_s+n_s*n_p);      % the sensitiity coordinates
vv      = reshape(v40,[n_p,n_s])';   % vector 2 matrix 
SenEq   = dfdx*vv+dfda;              % SENSITIVITY EQUATION
dy_sens = reshape(SenEq',n_s*n_p,1); % matrix 2 vector

dy      = [dy_state; dy_sens];            % returning states + sens

end

function dy = rhs_sens_methane(t,x,pars)

DG0     = pars(1);
R       = pars(2);
k       = pars(3);
nup     = pars(4);
DGp     = pars(5);
chi     = pars(6);
Y       = pars(7);
T       = pars(8);
Kac     = pars(9);
%Kn      = pars(10);
m       = pars(10);


   % specify auxiliary equations
    r = k*x(4)*x(1)/(Kac + x(1))*...
        (1 - exp(...
        (DG0+R*T*log(x(2)*x(3)/x(1))+nup*DGp)/...
        chi/R/T));
    
 %r = @(x) p.k*x(4)*x(1)/(p.Kac + x(1));
 
 % specify system of ODEs
    dy=zeros(4,1);
    
    dy(1) =-r;
    dy(2) = r;
    dy(3) = r;
    dy(4) = Y*r - m*x(4);  
    
% sensitivities    
    
% s11=x(5);
% s12=x(6);

% s21=x(7);
% s22=x(8);

% s31=x(9);
% s32=x(10);

% s41=x(11);
% s42=x(12);

s11=x(5);
s12=x(6);
s13=x(7);
s14=x(8);
s15=x(9);
s16=x(10);

s21=x(11);
s22=x(12);
s23=x(13);
s24=x(14);
s25=x(15);
s26=x(16);

s31=x(17);
s32=x(18);
s33=x(19);
s34=x(20);
s35=x(21);
s36=x(22);

s41=x(23);
s42=x(24);
s43=x(25);
s44=x(26);
s45=x(27);
s46=x(28);

xb=x(1);
y =x(2);
z =x(3);
u =x(4);


fexp=exp((DG0 + DGp*nup)/(chi*R*T));
fchi=((y*z)/xb)^(1/chi);
flog = log(y*z/xb);


dy(5)=(chi*((-1 + fchi*fexp)*Kac*(s11*u + s41*xb) + xb*(u - fchi*fexp*u + (-1 + fchi*fexp)*s41*xb))*y*z + fchi*fexp*u*(Kac + xb)*(s31*xb*y + s21*xb*z - s11*y*z))/(chi*(Kac + xb)^2*y*z);

dy(6)=(2*(chi*(-1 + fchi*fexp)*(Kac*s12*u + Kac*s42*xb + s42*xb^2)*y*z + fchi*fexp*u*(Kac + xb)*(s32*xb*y + s22*xb*z - s12*y*z)))/(chi*(Kac + xb)^2*y*z);

dy(7)=(3*(chi*(-1 + fchi*fexp)*(Kac*s13*u + Kac*s43*xb + s43*xb^2)*y*z + fchi*fexp*u*(Kac + xb)*(s33*xb*y + s23*xb*z - s13*y*z)))/(chi*(Kac + xb)^2*y*z);

dy(8)=(4*(chi*(-1 + fchi*fexp)*R*T*(Kac*s14*u + Kac*s44*xb + s44*xb^2)*y*z + fchi*u*(Kac + xb)*(DGp*fexp*xb*y*z + fexp*R*T*(s34*xb*y + s24*xb*z - s14*y*z))))/(chi*R*T*(Kac + xb)^2*y*z);

dy(9)=(chi*((-1 + fchi*fexp)*u*xb*(Kac + xb) + 5*(-1 + fchi*fexp)*(Kac*s15*u + Kac*s45*xb + s45*xb^2))*y*z + 5*fchi*fexp*u*(Kac + xb)*(s35*xb*y + s25*xb*z - s15*y*z))/(chi*(Kac + xb)^2*y*z);

dy(10)=(6*((s16*u*(-(chi*Kac) + fchi*fexp*((-1 + chi)*Kac - xb)))/chi + (-1 + fchi*fexp)*s46*xb*(Kac + xb) - (fchi*fexp*(DG0 + DGp*nup + flog*R*T)*u*xb*(Kac + xb))/(chi^2*R*T) + (fexp*s36*u*(Kac + xb)*y*((y*z)/xb)^(-1 + chi^(-1)))/chi + (fexp*s26*u*(Kac + xb)*z*((y*z)/xb)^(-1 + chi^(-1)))/chi))/(Kac + xb)^2;

dy(11)=-((chi*(-1 + fchi*fexp)*(Kac*s11*u + Kac*s41*xb - u*xb + s41*xb^2)*y*z + fchi*fexp*u*(Kac + xb)*(s31*xb*y + s21*xb*z - s11*y*z))/(chi*(Kac + xb)^2*y*z));

dy(12)=(-2*(chi*(-1 + fchi*fexp)*(Kac*s12*u + Kac*s42*xb + s42*xb^2)*y*z + fchi*fexp*u*(Kac + xb)*(s32*xb*y + s22*xb*z - s12*y*z)))/(chi*(Kac + xb)^2*y*z);

dy(13)=(-3*(chi*(-1 + fchi*fexp)*(Kac*s13*u + Kac*s43*xb + s43*xb^2)*y*z + fchi*fexp*u*(Kac + xb)*(s33*xb*y + s23*xb*z - s13*y*z)))/(chi*(Kac + xb)^2*y*z);

dy(14)=(-4*(chi*(-1 + fchi*fexp)*R*T*(Kac*s14*u + Kac*s44*xb + s44*xb^2)*y*z + fchi*fexp*u*(Kac + xb)*(DGp*xb*y*z + R*T*(s34*xb*y + s24*xb*z - s14*y*z))))/(chi*R*T*(Kac + xb)^2*y*z);

dy(15)=-((chi*(-1 + fchi*fexp)*(u*xb*(Kac + xb) + 5*(Kac*s15*u + Kac*s45*xb + s45*xb^2))*y*z + 5*fchi*fexp*u*(Kac + xb)*(s35*xb*y + s25*xb*z - s15*y*z))/(chi*(Kac + xb)^2*y*z));

dy(16)=(6*((1 - fchi*fexp)*s46*xb*(Kac + xb) + (fchi*fexp*(DG0 + DGp*nup + flog*R*T)*u*xb*(Kac + xb))/(chi^2*R*T) + (s16*u*(chi*Kac + fchi*fexp*(Kac - chi*Kac + xb)))/chi - (fexp*s36*u*(Kac + xb)*y*((y*z)/xb)^(-1 + chi^(-1)))/chi - (fexp*s26*u*(Kac + xb)*z*((y*z)/xb)^(-1 + chi^(-1)))/chi))/(Kac + xb)^2;

dy(17)=-((chi*(-1 + fchi*fexp)*(Kac*s11*u + Kac*s41*xb - u*xb + s41*xb^2)*y*z + fchi*fexp*u*(Kac + xb)*(s31*xb*y + s21*xb*z - s11*y*z))/(chi*(Kac + xb)^2*y*z));

dy(18)=(-2*(chi*(-1 + fchi*fexp)*(Kac*s12*u + Kac*s42*xb + s42*xb^2)*y*z + fchi*fexp*u*(Kac + xb)*(s32*xb*y + s22*xb*z - s12*y*z)))/(chi*(Kac + xb)^2*y*z);

dy(19)=(-3*(chi*(-1 + fchi*fexp)*(Kac*s13*u + Kac*s43*xb + s43*xb^2)*y*z + fchi*fexp*u*(Kac + xb)*(s33*xb*y + s23*xb*z - s13*y*z)))/(chi*(Kac + xb)^2*y*z);

dy(20)=(-4*(chi*(-1 + fchi*fexp)*R*T*(Kac*s14*u + Kac*s44*xb + s44*xb^2)*y*z + fchi*fexp*u*(Kac + xb)*(DGp*xb*y*z + R*T*(s34*xb*y + s24*xb*z - s14*y*z))))/(chi*R*T*(Kac + xb)^2*y*z);

dy(21)=-((chi*(-1 + fchi*fexp)*(u*xb*(Kac + xb) + 5*(Kac*s15*u + Kac*s45*xb + s45*xb^2))*y*z + 5*fchi*fexp*u*(Kac + xb)*(s35*xb*y + s25*xb*z - s15*y*z))/(chi*(Kac + xb)^2*y*z));

dy(22)=(6*((1 - fchi*fexp)*s46*xb*(Kac + xb) + (fchi*fexp*(DG0 + DGp*nup + flog*R*T)*u*xb*(Kac + xb))/(chi^2*R*T) + (s16*u*(chi*Kac + fchi*fexp*(Kac - chi*Kac + xb)))/chi - (fexp*s36*u*(Kac + xb)*y*((y*z)/xb)^(-1 + chi^(-1)))/chi - (fexp*s26*u*(Kac + xb)*z*((y*z)/xb)^(-1 + chi^(-1)))/chi))/(Kac + xb)^2;

dy(23)=-((chi*y*(Kac^2*m*s41 + xb*(m*s41*xb - (-1 + fchi*fexp)*(u - s41*xb)*Y) + Kac*(2*m*s41*xb + (-1 + fchi*fexp)*(s11*u + s41*xb)*Y))*z + fchi*fexp*u*(Kac + xb)*Y*(s31*xb*y + s21*xb*z - s11*y*z))/(chi*(Kac + xb)^2*y*z));

dy(24)=-((s42*(Kac + xb)*(Kac*m + xb*(m + 2*(-1 + fchi*fexp)*Y)) + (u*(chi*y*(Kac^2 + xb^2 + Kac*(2*xb + 2*(-1 + fchi*fexp)*s12*Y))*z + 2*fchi*fexp*(Kac + xb)*Y*(s32*xb*y + s22*xb*z - s12*y*z)))/(chi*y*z))/(Kac + xb)^2);

dy(25)=-((chi*y*(Kac^2*m*s43 + xb^2*(m*s43 + 3*(-1 + fchi*fexp)*(u + s43*Y)) + Kac*(2*m*s43*xb + 3*(-1 + fchi*fexp)*(s43*xb*Y + u*(xb + s13*Y))))*z + 3*fchi*fexp*u*(Kac + xb)*Y*(s33*xb*y + s23*xb*z - s13*y*z))/(chi*(Kac + xb)^2*y*z));

dy(26)=-(((4*DGp*fchi*fexp*u*xb*(Kac + xb)*Y)/(chi*R*T) - (4*s14*u*(chi*Kac + fchi*fexp*(Kac - chi*Kac + xb))*Y)/chi + s44*(Kac + xb)*(Kac*m + xb*(m + 4*(-1 + fchi*fexp)*Y)) + (4*fexp*s34*u*(Kac + xb)*y*Y*((y*z)/xb)^(-1 + chi^(-1)))/chi + (4*fexp*s24*u*(Kac + xb)*Y*z*((y*z)/xb)^(-1 + chi^(-1)))/chi)/(Kac + xb)^2);

dy(27)=-((chi*y*(Kac^2*m*s45 + xb^2*(m*s45 + (-1 + fchi*fexp)*(5*s45 + u)*Y) + Kac*(2*m*s45*xb + (-1 + fchi*fexp)*(5*s15*u + 5*s45*xb + u*xb)*Y))*z + 5*fchi*fexp*u*(Kac + xb)*Y*(s35*xb*y + s25*xb*z - s15*y*z))/(chi*(Kac + xb)^2*y*z));

dy(28)=((6*fchi*fexp*(DG0 + DGp*nup + flog*R*T)*u*xb*(Kac + xb)*Y)/(chi^2*R*T) + (6*s16*u*(chi*Kac + fchi*fexp*(Kac - chi*Kac + xb))*Y)/chi - s46*(Kac + xb)*(Kac*m + xb*(m + 6*(-1 + fchi*fexp)*Y)) - (6*fexp*s36*u*(Kac + xb)*y*Y*((y*z)/xb)^(-1 + chi^(-1)))/chi - (6*fexp*s26*u*(Kac + xb)*Y*z*((y*z)/xb)^(-1 + chi^(-1)))/chi)/(Kac + xb)^2;
end
