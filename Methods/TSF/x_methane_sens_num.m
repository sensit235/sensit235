function x_methane_sens_num

close all

addpath DERIVESTsuite/

% DATA
t_data=[0 2 4 5 6 10 12 14 16 18]*24*60*60;
y_data=[3.68 2.54 1.90 1.85 1.61 1.31 1.28 1.17 1.02 0.97]*1e-3;

% SIMULATION INTERVAL
days = 18;
td=(0:0.01:1)*days*24*60*60;
n_s = 4;    % number of states
n_p = 10;   % number of parameters

% PARAMETERS
DG0 = -15802.1961;  % standard value of Gibbs free energy at T=310.15K (J/mol)
R   = 8.3145;       % gas constant (J/mol/K)
k   = 2.5e-6;       % (mol/g/s) **
nup = 0.5;          % (per reaction) **
DGp = 45000;        % phosphorylation energy (J/(mol ATP)) **
chi = 2;            % (per reaction) **
Y   = 2.1;          % (Ymax) (g/mol) **
T   = 310.15;       % physiological temperature (K)
Kac = 5e-3;         % from Chewy  % (Kd) half saturation constant (molal) **
m   = 2.2e-7;       % from Chewy      % (D) specific maintenance rate (1/s) **
%Kn  = 0;           % effect of nutrient (molal)
parms = [DG0 R k nup DGp chi Y T Kac m];


% T_SPAN
tspan = [td(1),td(end)];

% INITIAL CONDITIONS
x0_fix      = [3.5e-3 45.2e-3 1e-6 0.009]; 

% OPTIONS
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

% SENS wrt PARS ANALYTICAL
x0_anal = [x0_fix zeros(1,24)];
[t,y]   = ode15s(@rhs_SenPar_anal,tspan,x0_anal,options,parms);

% SENS wrt PARs NUMERICAL
x0_num = [x0_fix zeros(1,n_s*n_p)];
[t_num,y_num] = ode15s(@rhs_SenPar_num,tspan,x0_num,options,parms,n_s,n_p);

% SENS wrt IC NUMERICAL
x0_num_IC = [x0_fix zeros(1,n_s*n_s)];
[t_ic,y_ic] = ode15s(@rhs_SenIC_num,tspan,x0_num_IC,options,parms,n_s);

% SAVE
save res_num.mat t y t_num y_num t_ic y_ic
end

% SENSITIVITIES WRT PARAMETERS
function dy = rhs_SenPar_num(t,v,pars,n_s,n_p)
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

% Model
dy_state = ss([v(1:n_s)],pars); 

% Jacobian df/dx
dfdx = jacobianest(@(x) ss(x,pars),[v(1:n_s)]);

% dF/da
dfda = jacobianest(@(x) ss([v(1:n_s)],x),pars);

% Sensitivity equation    
v40     = v(n_s+1:n_s+n_s*n_p);      % the sensitiity coordinates
vv      = reshape(v40,[n_p,n_s])';   % vector 2 matrix 
SenEq   = dfdx*vv+dfda;              % SENSITIVITY EQUATION
dy_sens = reshape(SenEq',n_s*n_p,1); % matrix 2 vector

dy      = [dy_state; dy_sens];            % returning states + sens

end

% SENSITIVITIES WRT INITIAL CONDITIONS
function dy = rhs_SenIC_num(t,v,pars,n_s)

rr = @(x,p) p(3)*x(4)*x(1)/(p(9) + x(1))*...
     (1 - exp(...
     (p(1)+p(2)*p(8)*log(x(2)*x(3)/x(1))+p(4)*p(5))/...
     p(6)/p(2)/p(8)));
 
ss = @(x,p)[-rr(x,p);... 
             rr(x,p);...
             rr(x,p);...
             p(7)*rr(x,p) - p(10)*x(4)];

% Model
dy_state = ss(v(1:n_s),pars); 

% Jacobian df/dx
dfdx = jacobianest(@(x) ss(x,pars),v(1:n_s));

% dF/da
%dfda = jacobianest(@(x) ss(v(1:n_s),x),pars);

% Sensitivity equation    
v_IC       = v(n_s+1:n_s+n_s*n_s);          % the sensitiity coordinates
vv_IC      = reshape(v_IC,[n_s,n_s])';      % vector 2 matrix 
Sen_IC     = dfdx*vv_IC;                    % SENSITIVITY EQUATION
dy_sens_IC = reshape(Sen_IC',n_s*n_s,1);    % matrix 2 vector
dy         = [dy_state; dy_sens_IC];        % returning states + sens_IC       
end


function dy = rhs_SenPar_anal(t,x,pars)

% pars = [DG0 R k nup DGp chi Y T Kac m];

DG0     = pars(1);
R       = pars(2);
k       = pars(3);
nup     = pars(4);
DGp     = pars(5);
chi     = pars(6);
Y       = pars(7);
T       = pars(8);
Kac     = pars(9);
m       = pars(10);


   % specify auxiliary equations
    r = k*x(4)*x(1)/(Kac + x(1))*(1 - exp((DG0+R*T*log(x(2)*x(3)/x(1))+nup*DGp)/chi/R/T));
    
 % specify system of ODEs
    dy=zeros(4,1);
    
    dy(1) =-r;
    dy(2) = r;
    dy(3) = r;
    dy(4) = Y*r - m*x(4);  
    
% STATE VARIABLES    
xb = x(1);
y  = x(2);
z  = x(3);
u  = x(4);
% SENSITIVITIES    
s11 = x(5);
s12 = x(6);
s13 = x(7);
s14 = x(8);
s15 = x(9);
s16 = x(10);

s21 = x(11);
s22 = x(12);
s23 = x(13);
s24 = x(14);
s25 = x(15);
s26 = x(16);

s31 = x(17);
s32 = x(18);
s33 = x(19);
s34 = x(20);
s35 = x(21);
s36 = x(22);

s41 = x(23);
s42 = x(24);
s43 = x(25);
s44 = x(26);
s45 = x(27);
s46 = x(28);

fexp=exp((DG0 + DGp*nup)/(chi*R*T));
fchi=((y*z)/xb)^(1/chi);
flog = log(y*z/xb);
fexp2=exp(-((-DG0 - DGp*nup - flog*R*T)/(chi*R*T)));

dy(5)=((1 - fexp2)*u*xb)/(Kac + xb)^2 - ((1 - fexp2)*s41*xb)/(Kac + xb) + s11*(((1 - fexp2)*u*xb)/(Kac + xb)^2 - ((1 - fexp2)*u)/(Kac + xb) - (fexp2*u)/(chi*(Kac + xb))) + (fexp2*s21*u*xb)/(chi*(Kac + xb)*y) + (fexp2*s31*u*xb)/(chi*(Kac + xb)*z);

dy(6)=(-2*(1 - fexp2)*s42*xb)/(Kac + xb) + s12*((2*(1 - fexp2)*u*xb)/(Kac + xb)^2 - (2*(1 - fexp2)*u)/(Kac + xb) - (2*fexp2*u)/(chi*(Kac + xb))) + (2*fexp2*s22*u*xb)/(chi*(Kac + xb)*y) + (2*fexp2*s32*u*xb)/(chi*(Kac + xb)*z);

dy(7)=(-3*(1 - fexp2)*s43*xb)/(Kac + xb) + s13*((3*(1 - fexp2)*u*xb)/(Kac + xb)^2 - (3*(1 - fexp2)*u)/(Kac + xb) - (3*fexp2*u)/(chi*(Kac + xb))) + (3*fexp2*s23*u*xb)/(chi*(Kac + xb)*y) + (3*fexp2*s33*u*xb)/(chi*(Kac + xb)*z);

dy(8)=(-4*(1 - fexp2)*s44*xb)/(Kac + xb) + (4*DGp*fexp2*u*xb)/(chi*R*T*(Kac + xb)) + s14*((4*(1 - fexp2)*u*xb)/(Kac + xb)^2 - (4*(1 - fexp2)*u)/(Kac + xb) - (4*fexp2*u)/(chi*(Kac + xb))) + (4*fexp2*s24*u*xb)/(chi*(Kac + xb)*y) + (4*fexp2*s34*u*xb)/(chi*(Kac + xb)*z);

dy(9)=(-5*(1 - fexp2)*s45*xb)/(Kac + xb) - ((1 - fexp2)*u*xb)/(Kac + xb) + s15*((5*(1 - fexp2)*u*xb)/(Kac + xb)^2 - (5*(1 - fexp2)*u)/(Kac + xb) - (5*fexp2*u)/(chi*(Kac + xb))) + (5*fexp2*s25*u*xb)/(chi*(Kac + xb)*y) + (5*fexp2*s35*u*xb)/(chi*(Kac + xb)*z);

dy(10)=(-6*(1 - fexp2)*s46*xb)/(Kac + xb) + (6*fexp2*(-DG0 - DGp*nup - flog*R*T)*u*xb)/(chi^2*R*T*(Kac + xb)) + s16*((6*(1 - fexp2)*u*xb)/(Kac + xb)^2 - (6*(1 - fexp2)*u)/(Kac + xb) - (6*fexp2*u)/(chi*(Kac + xb))) + (6*fexp2*s26*u*xb)/(chi*(Kac + xb)*y) + (6*fexp2*s36*u*xb)/(chi*(Kac + xb)*z);

dy(11)=-(((1 - fexp2)*u*xb)/(Kac + xb)^2) + ((1 - fexp2)*s41*xb)/(Kac + xb) + s11*(-(((1 - fexp2)*u*xb)/(Kac + xb)^2) + ((1 - fexp2)*u)/(Kac + xb) + (fexp2*u)/(chi*(Kac + xb))) - (fexp2*s21*u*xb)/(chi*(Kac + xb)*y) - (fexp2*s31*u*xb)/(chi*(Kac + xb)*z);

dy(12)=(2*(1 - fexp2)*s42*xb)/(Kac + xb) + s12*((-2*(1 - fexp2)*u*xb)/(Kac + xb)^2 + (2*(1 - fexp2)*u)/(Kac + xb) + (2*fexp2*u)/(chi*(Kac + xb))) - (2*fexp2*s22*u*xb)/(chi*(Kac + xb)*y) - (2*fexp2*s32*u*xb)/(chi*(Kac + xb)*z);

dy(13)=(3*(1 - fexp2)*s43*xb)/(Kac + xb) + s13*((-3*(1 - fexp2)*u*xb)/(Kac + xb)^2 + (3*(1 - fexp2)*u)/(Kac + xb) + (3*fexp2*u)/(chi*(Kac + xb))) - (3*fexp2*s23*u*xb)/(chi*(Kac + xb)*y) - (3*fexp2*s33*u*xb)/(chi*(Kac + xb)*z);

dy(14)=(4*(1 - fexp2)*s44*xb)/(Kac + xb) - (4*DGp*fexp2*u*xb)/(chi*R*T*(Kac + xb)) + s14*((-4*(1 - fexp2)*u*xb)/(Kac + xb)^2 + (4*(1 - fexp2)*u)/(Kac + xb) + (4*fexp2*u)/(chi*(Kac + xb))) - (4*fexp2*s24*u*xb)/(chi*(Kac + xb)*y) - (4*fexp2*s34*u*xb)/(chi*(Kac + xb)*z);

dy(15)=(5*(1 - fexp2)*s45*xb)/(Kac + xb) + ((1 - fexp2)*u*xb)/(Kac + xb) + s15*((-5*(1 - fexp2)*u*xb)/(Kac + xb)^2 + (5*(1 - fexp2)*u)/(Kac + xb) + (5*fexp2*u)/(chi*(Kac + xb))) - (5*fexp2*s25*u*xb)/(chi*(Kac + xb)*y) - (5*fexp2*s35*u*xb)/(chi*(Kac + xb)*z);

dy(16)=(6*(1 - fexp2)*s46*xb)/(Kac + xb) - (6*fexp2*(-DG0 - DGp*nup - flog*R*T)*u*xb)/(chi^2*R*T*(Kac + xb)) + s16*((-6*(1 - fexp2)*u*xb)/(Kac + xb)^2 + (6*(1 - fexp2)*u)/(Kac + xb) + (6*fexp2*u)/(chi*(Kac + xb))) - (6*fexp2*s26*u*xb)/(chi*(Kac + xb)*y) - (6*fexp2*s36*u*xb)/(chi*(Kac + xb)*z);

dy(17)=-(((1 - fexp2)*u*xb)/(Kac + xb)^2) + ((1 - fexp2)*s41*xb)/(Kac + xb) + s11*(-(((1 - fexp2)*u*xb)/(Kac + xb)^2) + ((1 - fexp2)*u)/(Kac + xb) + (fexp2*u)/(chi*(Kac + xb))) - (fexp2*s21*u*xb)/(chi*(Kac + xb)*y) - (fexp2*s31*u*xb)/(chi*(Kac + xb)*z);

dy(18)=(2*(1 - fexp2)*s42*xb)/(Kac + xb) + s12*((-2*(1 - fexp2)*u*xb)/(Kac + xb)^2 + (2*(1 - fexp2)*u)/(Kac + xb) + (2*fexp2*u)/(chi*(Kac + xb))) - (2*fexp2*s22*u*xb)/(chi*(Kac + xb)*y) - (2*fexp2*s32*u*xb)/(chi*(Kac + xb)*z);

dy(19)=(3*(1 - fexp2)*s43*xb)/(Kac + xb) + s13*((-3*(1 - fexp2)*u*xb)/(Kac + xb)^2 + (3*(1 - fexp2)*u)/(Kac + xb) + (3*fexp2*u)/(chi*(Kac + xb))) - (3*fexp2*s23*u*xb)/(chi*(Kac + xb)*y) - (3*fexp2*s33*u*xb)/(chi*(Kac + xb)*z);

dy(20)=(4*(1 - fexp2)*s44*xb)/(Kac + xb) - (4*DGp*fexp2*u*xb)/(chi*R*T*(Kac + xb)) + s14*((-4*(1 - fexp2)*u*xb)/(Kac + xb)^2 + (4*(1 - fexp2)*u)/(Kac + xb) + (4*fexp2*u)/(chi*(Kac + xb))) - (4*fexp2*s24*u*xb)/(chi*(Kac + xb)*y) - (4*fexp2*s34*u*xb)/(chi*(Kac + xb)*z);

dy(21)=(5*(1 - fexp2)*s45*xb)/(Kac + xb) + ((1 - fexp2)*u*xb)/(Kac + xb) + s15*((-5*(1 - fexp2)*u*xb)/(Kac + xb)^2 + (5*(1 - fexp2)*u)/(Kac + xb) + (5*fexp2*u)/(chi*(Kac + xb))) - (5*fexp2*s25*u*xb)/(chi*(Kac + xb)*y) - (5*fexp2*s35*u*xb)/(chi*(Kac + xb)*z);

dy(22)=(6*(1 - fexp2)*s46*xb)/(Kac + xb) - (6*fexp2*(-DG0 - DGp*nup - flog*R*T)*u*xb)/(chi^2*R*T*(Kac + xb)) + s16*((-6*(1 - fexp2)*u*xb)/(Kac + xb)^2 + (6*(1 - fexp2)*u)/(Kac + xb) + (6*fexp2*u)/(chi*(Kac + xb))) - (6*fexp2*s26*u*xb)/(chi*(Kac + xb)*y) - (6*fexp2*s36*u*xb)/(chi*(Kac + xb)*z);

dy(23)=-(((1 - fexp2)*u*xb*Y)/(Kac + xb)^2) - (fexp2*s21*u*xb*Y)/(chi*(Kac + xb)*y) + s11*(-(((1 - fexp2)*u*xb*Y)/(Kac + xb)^2) + ((1 - fexp2)*u*Y)/(Kac + xb) + (fexp2*u*Y)/(chi*(Kac + xb))) + s41*(-m + ((1 - fexp2)*xb*Y)/(Kac + xb)) - (fexp2*s31*u*xb*Y)/(chi*(Kac + xb)*z);

dy(24)=-u - (2*fexp2*s22*u*xb*Y)/(chi*(Kac + xb)*y) + s12*((-2*(1 - fexp2)*u*xb*Y)/(Kac + xb)^2 + (2*(1 - fexp2)*u*Y)/(Kac + xb) + (2*fexp2*u*Y)/(chi*(Kac + xb))) + s42*(-m + (2*(1 - fexp2)*xb*Y)/(Kac + xb)) - (2*fexp2*s32*u*xb*Y)/(chi*(Kac + xb)*z);

dy(25)=(3*(1 - fexp2)*u*xb)/(Kac + xb) - (3*fexp2*s23*u*xb*Y)/(chi*(Kac + xb)*y) + s13*((-3*(1 - fexp2)*u*xb*Y)/(Kac + xb)^2 + (3*(1 - fexp2)*u*Y)/(Kac + xb) + (3*fexp2*u*Y)/(chi*(Kac + xb))) + s43*(-m + (3*(1 - fexp2)*xb*Y)/(Kac + xb)) - (3*fexp2*s33*u*xb*Y)/(chi*(Kac + xb)*z);

dy(26)=(-4*DGp*fexp2*u*xb*Y)/(chi*R*T*(Kac + xb)) - (4*fexp2*s24*u*xb*Y)/(chi*(Kac + xb)*y) + s14*((-4*(1 - fexp2)*u*xb*Y)/(Kac + xb)^2 + (4*(1 - fexp2)*u*Y)/(Kac + xb) + (4*fexp2*u*Y)/(chi*(Kac + xb))) + s44*(-m + (4*(1 - fexp2)*xb*Y)/(Kac + xb)) - (4*fexp2*s34*u*xb*Y)/(chi*(Kac + xb)*z);

dy(27)=((1 - fexp2)*u*xb*Y)/(Kac + xb) - (5*fexp2*s25*u*xb*Y)/(chi*(Kac + xb)*y) + s15*((-5*(1 - fexp2)*u*xb*Y)/(Kac + xb)^2 + (5*(1 - fexp2)*u*Y)/(Kac + xb) + (5*fexp2*u*Y)/(chi*(Kac + xb))) + s45*(-m + (5*(1 - fexp2)*xb*Y)/(Kac + xb)) - (5*fexp2*s35*u*xb*Y)/(chi*(Kac + xb)*z);

dy(28)=(-6*fexp2*(-DG0 - DGp*nup - flog*R*T)*u*xb*Y)/(chi^2*R*T*(Kac + xb)) - (6*fexp2*s26*u*xb*Y)/(chi*(Kac + xb)*y) + s16*((-6*(1 - fexp2)*u*xb*Y)/(Kac + xb)^2 + (6*(1 - fexp2)*u*Y)/(Kac + xb) + (6*fexp2*u*Y)/(chi*(Kac + xb))) + s46*(-m + (6*(1 - fexp2)*xb*Y)/(Kac + xb)) - (6*fexp2*s36*u*xb*Y)/(chi*(Kac + xb)*z);



end
