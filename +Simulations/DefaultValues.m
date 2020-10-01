% CREDITS
% Made by Rahul S. Yerrabelli under Nai-Kong Cheung, MD, PhD at
% Memorial Sloan Kettering Cancer Center
% This is part of the code accompanying the publication of "Intra-Ommaya
% Compartmental Radioimmunotherapy using 131I-Omburtamab? Pharmacokinetic
% Modeling to Optimize Therapeutic Index" in European Journal of Nuclear
% Medicine and Molecular Imaging (EJNMMI) (Accepted Sep 19, 2020). See that
% paper and its supplement file for a greater description of the model and
% results.


% DESCRIPTION
% This file sets the default values used for the optimization. This file
% should be run before every optimization to ensure the default values are
% used. Of course, each optimization procedure can have specific values
% that are different than these values (i.e. to see how changing the value
% from default would affect results), but that should be done in the
% optimization code after running this code.


% DETAILS
% Tested in MATLAB® 2017a for Mac.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The amount of activity administered
activityAdmin = 50; %mCi (note: 1mCi=37MBq)


% Fitted anatomical values, taken from the median of allfits lsq0 skipping 
% increasing beginning
Vv = 36.611/1000; %Vv is ultimately in L (thus /1000 for mL->L)
CL = 11.972/1000/3600; %CL is ultimately in L/s (thus /1000/3600 for mL/hr->K/s)


% time_values_needed is ultimately in sec
% time_values_needed = [1 3 12 24 48]*3600; %up to 48 hours
time_values_needed = [0:1:1440*100]*60; %100 days worth


% Fraction of antibody that is still immunoreactive (and thus useful) after
% radiolabeling
immunoreactivity = 0.69; %doesn't really matter for fitting part


N_A=6.02214085774e23; % Avogadro's number. last two digits unsure


% tumorload is the initial [receptor] in mol/L. R0 is the same, but in
% antigens/mL
% He et al's R0 for 3F8/GD2 was 2.8635e+14 ant/mL -> 4.755e-7 mol/L
% I calculated 8H9/B7-H3's R0 assuming fcov=3%, fexp=50%, D_T=10um, and the
% NR value (antigen/cell) from the median of the sampled data=41667.1162040375
% As a normal value, use:
tumorload = 1.193674951486256e+12 / (6.022e23) *1000; %mol/L
% Using other NR's besides the median, I got a 10-90% data range of
% 3.9673e9 to 3.179e13; complete range: 2.4165e9 to 6.4754e13 ant/mL.

% If you wanted to actually calculate the R0 in code, it would be
%fexp=0.5;
%fcov=0.3;
%D_T=10e-4; %e-4 to go from um to cm
%C_SC=1/(4/3*pi*(D_T/2)^3);
%N_R=41667;
%tumorload = fexp*fcov*C_SC*N_R  * 1/N_A*1000;


perTV = 0; %percent of tumor cells that are in the ventricles
perNB = 0.01; %amount on nonspecific binding (%) relative to speecific binding
kAR = 1.846e4; %+/-.006e4, which is specifically for m8H9 with 4Ig-B7-H3
k_AR = 1.650e-4;%+/-.002e-4, which is specifically for m8H9 with 4Ig-B7-H3
V = 0.140; %total volume, doesn't really matter for fitting part (liters)
Vs = V - Vv; %Vs = Subarachnoid space volume (liters)
n = 1; %Number of layers of cells lining the cavity. Never tested values outside of 1. (unitless) 
MM = 147507; %Molar mass, Da (2.4908102999999997e-19 * 6.022e23)
cI0 = 50; %specific activity, mCi/mg or Ci/g (50mCi/mg=185MBq/mg)
%t_half = half life
t_half = 193*3600; %convert hours to seconds


% Surface Area of entire CSF space
% Scaled knowing cited value that 140mL CSF volume has 0.18m2 (1800cm2)
% Multiply by 100 to convert below units from m2 to dm2 so that units
% become L when you multiply with D_T
S = 100*0.18*sqrt((Vv+Vs)./0.14); %units are m2

