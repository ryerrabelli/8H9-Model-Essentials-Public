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
% This function represents the two compartment 8H9 pharmacokinetic model.
% It contains the 12 autonomous ordinary differential equations (ODE) that 
% define it. It can be solved by an ODE solver in MATLAB® (i.e. ode45).
% MATLAB views a system of multiple ODEs as one ODE with vector inputs and
% vector outputs. Thus the input and output are each a 12-length vector.


% ARGUMENTS
% t:        Time. As ODEs are autonomous, not relevant for calculation.
% y:        Vector of ODE inputs. (Units: M=mol/L)
% immunoreactivity:      
%           Immunoreactivity (unitless, on [0,1))
% perTV:    TumorVent, percent of tumor in ventricles (unitless, on [0,1])
% perNB     Percent nonspecificity (unitless, [0,1])
% kAR:      Association rate constant of forming antigen-receptor complex
%           (units: 1/(M s))
% k_AR:     Dissociation rate constant of forming antigen-receptor complex
%           (units: 1/s)
% Vv:       Ventricular volume (units: liters)
% Vs:       Subarachnoid space volume (units: liters)
% n         Number of layers of cells lining the cavity. Never tested 
%           values outside of 1. (unitless) 
% CL        Clearance (L/s)
% S_dm2     Surface area of entire CSF space (units: dm^2)
% D_T       Cell diameter (units: dm)


% OUTPUTS
% dy:       A vector of outputs where each element is the output of one of
%           the autonomous ODEs at the given input y. (Units: M/s)


% DETAILS
% Tested in MATLAB® 2017a for Mac.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dy = ModOde2(t, y,...
    immunoreactivity, perTV, perNB, kAR, k_AR, Vv, Vs, n, CL, S_dm2, D_T)


dy = zeros(12,1); %initialize output vector

% Split input concentration vector into labeled elements. Units in M=mol/L
cBV = y(1);
cNV = y(2);
cBS = y(3);
cNS = y(4);
cARV = y(5);
cRV = y(6);
cARS = y(7);
cRS = y(8);
cNARV = y(9);
cNRV = y(10);
cNARS = y(11);
cNRS = y(12);


% Parameters
k_NAR = k_AR;
kNAR = kAR * perNB;
INF = 0; %continuous infusion rate


% Differential Equations
dcBV = (INF * immunoreactivity / Vv) - (CL*cBV / Vv)...
    + perNB*(n*D_T*S_dm2/Vv)*(k_NAR*cNARV - kNAR*cBV*cNRV)...
    + perTV*(n*D_T*S_dm2/Vv)*(k_AR*cARV - kAR*cBV*cRV);
dcNV = (INF * (1 - immunoreactivity) / Vv) - (CL*cNV / Vv);

dcBS = (CL*cBV / Vs) - (CL*cBS / Vs)...
    + (1-perTV)*(n*D_T*S_dm2/Vs)*(k_AR*cARS - kAR*cBS*cRS)...
    + perNB*(n*D_T*S_dm2/Vs)*(k_NAR*cNARS - kNAR*cBS*cNRS);
dcNS = (CL*cNV / Vs) - (CL*cNS / Vs);

dcARV = kAR*cBV*cRV - k_AR*cARV;
dcRV = k_AR*cARV - kAR*cBV*cRV;

dcARS = kAR*cBS*cRS - k_AR*cARS;
dcRS = k_AR*cARS - kAR*cBS*cRS;

dcNARV = kNAR*cBV*cNRV - k_NAR*cNARV;
dcNRV = k_NAR*cNARV - kNAR*cBV*cNRV; %***** cRV? cNRV!

dcNARS = kNAR*cBS*cNRS - k_NAR*cNARS; %***cRS? cNRS!
dcNRS = k_NAR*cNARS - kNAR*cBS*cNRS;


% Export individual concentration rate variables into a vector of them
dy(1) = dcBV;
dy(2) = dcNV;
dy(3) = dcBS;
dy(4) = dcNS;
dy(5) = dcARV;
dy(6) = dcRV;
dy(7) = dcARS;
dy(8) = dcRS;
dy(9) = dcNARV;
dy(10) = dcNRV;
dy(11) = dcNARS;
dy(12) = dcNRS;