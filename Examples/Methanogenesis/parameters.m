%% parameters
%
% This script initialises a data structure with all of the parameters

p.DG0 = -15802.1961; % standard value of Gibbs free energy at T=310.15K (J/mol)
p.R = 8.3145; % gas constant (J/mol/K)
p.k = 2.5e-6; % (mol/g/s) **
p.nup = 0.5; % (per reaction) **
p.DGp = 45000; % phosphorylation energy (J/(mol ATP))
p.chi = 2; % (per reaction) **
p.Y = 2.1; % (Ymax) (g/mol) **
p.T = 310.15; % physiological temperature (K)
p.Kac = 5e-3; % (Kd) half saturation constant (molal) **
p.Kn = 0; % effect of nutrient (molal)
p.m = 2.2e-7; % (D) specific maintenance rate (1/s) **
p.dummy = 0.;

% minimum

% pmin.DG0 = -15802.1961; % standard value of Gibbs free energy at T=310.15K (J/mol)
% pmin.R = 8.3145; % gas constant (J/mol/K)
% pmin.k = 0.2e-6; % (mol/g/s) **
% pmin.nup = 0; % (per reaction) **
% pmin.DGp = 45000; % phosphorylation energy (J/(mol ATP))
% pmin.chi = 1; % (per reaction) **
% pmin.Y = 2.1; % (Ymax) (g/mol) **
% pmin.T = 310.15; % physiological temperature (K)
% pmin.Kac = 1e-3; % (Kd) half saturation constant (molal) **
% pmin.Kn = 0; % effect of nutrient (molal)
% pmin.m = 5.1e-8; % (D) specific maintenance rate (1/s) **

% pmin.DG0 = -15802.1961; % standard value of Gibbs free energy at T=310.15K (J/mol)
% pmin.R = 8.3145; % gas constant (J/mol/K)
% pmin.k = 0.2e-6; % (mol/g/s) **
% pmin.nup = 0.3; % (per reaction) **
% pmin.DGp = 45000; % phosphorylation energy (J/(mol ATP))
% pmin.chi = 1; % (per reaction) **
% pmin.Y = 2.1; % (Ymax) (g/mol) **
% pmin.T = 310.15; % physiological temperature (K)
% pmin.Kac = 1e-3; % (Kd) half saturation constant (molal) **
% pmin.Kn = 0; % effect of nutrient (molal)
% pmin.m = 5.1e-8; % (D) specific maintenance rate (1/s) **

pmin.DG0 = -15802.1961; % standard value of Gibbs free energy at T=310.15K (J/mol)
pmin.R = 8.3145; % gas constant (J/mol/K)
pmin.k = 0.2e-6; % (mol/g/s) **
pmin.nup = 0.25; % (per reaction) **
pmin.DGp = 45000; % phosphorylation energy (J/(mol ATP))
pmin.chi = 1; % (per reaction) **
pmin.Y = 2.1; % (Ymax) (g/mol) **
pmin.T = 310.15; % physiological temperature (K)
pmin.Kac = 1e-3; % (Kd) half saturation constant (molal) **
pmin.Kn = 0; % effect of nutrient (molal)
pmin.m = 1e-7; % (D) specific maintenance rate (1/s) **
pmin.dummy = 0.;

% maximum

% pmax.DG0 = -15802.1961; % standard value of Gibbs free energy at T=310.15K (J/mol)
% pmax.R = 8.3145; % gas constant (J/mol/K)
% pmax.k = 6.5e-6; % (mol/g/s) **
% pmax.nup = 1.0; % (per reaction) **
% pmax.DGp = 45000; % phosphorylation energy (J/(mol ATP))
% pmax.chi = 2; % (per reaction) **
% pmax.Y = 4; % (Ymax) (g/mol) **
% pmax.T = 310.15; % physiological temperature (K)
% pmax.Kac = 29.1e-3; % (Kd) half saturation constant (molal) **
% pmax.Kn = 0; % effect of nutrient (molal)
% pmax.m = 9.4e-5; % (D) specific maintenance rate (1/s) **

% pmax.DG0 = -15802.1961; % standard value of Gibbs free energy at T=310.15K (J/mol)
% pmax.R = 8.3145; % gas constant (J/mol/K)
% pmax.k = 6.5e-6; % (mol/g/s) **
% pmax.nup = 0.7; % (per reaction) **
% pmax.DGp = 45000; % phosphorylation energy (J/(mol ATP))
% pmax.chi = 2; % (per reaction) **
% pmax.Y = 4; % (Ymax) (g/mol) **
% pmax.T = 310.15; % physiological temperature (K)
% pmax.Kac = 29.1e-3; % (Kd) half saturation constant (molal) **
% pmax.Kn = 0; % effect of nutrient (molal)
% pmax.m = 9.4e-5; % (D) specific maintenance rate (1/s) **

pmax.DG0 = -15802.1961; % standard value of Gibbs free energy at T=310.15K (J/mol)
pmax.R = 8.3145; % gas constant (J/mol/K)
pmax.k = 6.5e-6; % (mol/g/s) **
pmax.nup = 0.75; % (per reaction) **
pmax.DGp = 45000; % phosphorylation energy (J/(mol ATP))
pmax.chi = 2; % (per reaction) **
pmax.Y = 4; % (Ymax) (g/mol) **
pmax.T = 310.15; % physiological temperature (K)
pmax.Kac = 29.1e-3; % (Kd) half saturation constant (molal) **
pmax.Kn = 0; % effect of nutrient (molal)
pmax.m = 3e-7; % (D) specific maintenance rate (1/s) **
pmax.dummy = 0.;