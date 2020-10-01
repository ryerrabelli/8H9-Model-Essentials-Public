%This function represents the two compartment 8H9 pharmacokinetic model

function dy = ModOde2(t, y,...
    eff, perTV, perNB, kAR, k_AR, Vv, Vs, n, CL, S, D_T)
%eff = immunoreactivity
%perTV = tumorvent
%perNB is specificity

dy = zeros(12,1);
%Concentrations, in mol/L
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

%Parameters
k_NAR = k_AR;
kNAR = kAR * perNB;
INF = 0; %continuous infusion rate


%Differential Equations
%perTV = TumorVent
dcBV = (INF * eff / Vv) - (CL*cBV / Vv)...
    + perNB*(n*D_T*S/Vv)*(k_NAR*cNARV - kNAR*cBV*cNRV)...
    + perTV*(n*D_T*S/Vv)*(k_AR*cARV - kAR*cBV*cRV);
dcNV = (INF * (1 - eff) / Vv) - (CL*cNV / Vv);

dcBS = (CL*cBV / Vs) - (CL*cBS / Vs)...
    + (1-perTV)*(n*D_T*S/Vs)*(k_AR*cARS - kAR*cBS*cRS)...
    + perNB*(n*D_T*S/Vs)*(k_NAR*cNARS - kNAR*cBS*cNRS);
dcNS = (CL*cNV / Vs) - (CL*cNS / Vs);

dcARV = kAR*cBV*cRV - k_AR*cARV;
dcRV = k_AR*cARV - kAR*cBV*cRV;

dcARS = kAR*cBS*cRS - k_AR*cARS;
dcRS = k_AR*cARS - kAR*cBS*cRS;

dcNARV = kNAR*cBV*cNRV - k_NAR*cNARV;
dcNRV = k_NAR*cNARV - kNAR*cBV*cNRV; %***** cRV? cNRV!

dcNARS = kNAR*cBS*cNRS - k_NAR*cNARS; %***cRS? cNRS!
dcNRS = k_NAR*cNARS - kNAR*cBS*cNRS;

%Rates
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