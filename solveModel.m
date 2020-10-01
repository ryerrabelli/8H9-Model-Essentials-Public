%**** IMPORTANT- If you change this, don't forget to change
%solveModelMultlInfus.m as well

% **** AUC changes based off of the time_values_needed
function [t, cIobs, AUC_cIAR, AUC_cIA, returns, AUC_rIAR, AUC_rIA, AUC_cIARV, AUC_cIARS]...
    = solveModel(activityAdmin, ...
    Vv, CL, time_values_needed, ... %time_values_needed is in seconds
    defined_params, return_all) 
%input S in m2 (will be converted to dm2 inside);

if (~exist('return_all', 'var'))
    return_all = false;
end

if (exist('defined_params', 'var'))
    immunoreactivity = defined_params(1);
    tumorload = defined_params(2); %aka R0
    perTV = defined_params(3); %TumorVent
    perNB = defined_params(4); %nonspecificity
    kAR = defined_params(5);
    k_AR = defined_params(6);
    V = defined_params(7); %L
    Vs = V - Vv; %L
    n = defined_params(8);
    MM = defined_params(9);
    cI0 = defined_params(10); %mCi/mg, Ci/g
    t_half = defined_params(11);
else
    disp('Warning: Using assumed values. V assumed to be 1000x normal');
    immunoreactivity = 0.69; %doesn't matter for fitting part
    tumorload = 4.755e-7; %TumorVent
    perTV = 0; %TumorVent, percent of tumor in ventricles%TumorVent, percent of tumor in ventricles
    perNB = 0.01; %nonspecificity
    kAR = 1.846e4; %+/-.006e4, which is specifically for m8H9 with 4Ig-B7-H3
    k_AR = 1.650e-4;%+/-.002e-4, which is specifically for m8H9 with 4Ig-B7-H3
    V = 0.140*1000; %doesn't matter for fitting part, make it 140,000mL
    Vs = V - Vv;
    n = 1;

    MM = 147507; %Da %2.4908102999999997e-19 * 6.022e23; %1499965.96266
    cI0 = 50; %mCi/mg, Ci/g, or 1850MBq/mg
    t_half = 193*3600; %convert hours to seconds

end

%&& has to be used, not & in order to use short-circuiting
if exist('defined_params', 'var') && length(defined_params) > 11
    S = defined_params(12) * 100; %convert from m2 to dm2
else
    %Surface Area of entire CSF space
    %Scaled knowing cited value that 140mL CSF volume has 0.18m2 (1800cm2)
    %Multiply by 100 to convert below units from m2 to dm2 so that units
    %become L when you multiply with D_T
    S = 100*0.18*sqrt((Vv+Vs)./0.14);
end

%units in dm so becomes L when multiply with surface area (in dm2)
D_T = 0.0001; %Tumor cell spherical diameter (0.0001dm = 10um)


%Initial Values
%B = free antibody w/ IR, N = free antibody w/o IR. See supplement.
%dose = (9.52e-8/2*0.140/Vv) * 3.84; %mol/L
cDose=activityAdmin/cI0*1/1000*1/MM*1/Vv; %activityAdmin in mCi => mol/L
cBV0 = immunoreactivity * cDose; %mol/L
cNV0 = (1-immunoreactivity) * cDose; %mol/L
cBS0 = 0;
cNS0 = 0;
cARV0 = 0;
%Starred one means "normal" is nonzero
cRV0 = perTV * tumorload; %mol/L,  tumorload is R0 (in moles/L) *
cARS0 = 0;
cRS0 = (1-perTV) * tumorload; %mol/L *
cNARV0 = 0;
cNRV0 = perNB * (Vv/V) * tumorload; %mol/L *
cNARS0 = 0;
cNRS0 = perNB * (Vs/V) * tumorload; %mol/L *

%Inputs
ini = [cBV0, cNV0, cBS0, cNS0, cARV0, cRV0, cARS0, cRS0, cNARV0, cNRV0, cNARS0, cNRS0];
%tspan = linspace(0,3000*60,60+1); %units is seconds with a point per minute
if time_values_needed(1) == 0
    tspan = time_values_needed;
elseif isnan(time_values_needed(1))
    %NaN start indicates ode isn't supposed to start from time 0 and thus
    %should bypass the check above
    tspan = [time_values_needed(2:end)];
else
    tspan = [0 time_values_needed];
end
%tolerance = 1e-16;
tolerance = 2.2205e-14;
if tspan(end) == Inf
    options = odeset('RelTol', tolerance, 'AbsTol', tolerance,...
        'Events', @(t,y)eventFcn(t,y,cDose)); %event stops the ode when values are low enough
else
    options = odeset('RelTol', tolerance, 'AbsTol', tolerance);
end


%*** NOW TIME TO SOLVE ODEs ***
%AbFinal(t, y, eff, perTV, perNB, kAR, k_AR, Vv, Vs, n, CL)
[t,y] = ode45(@(t,y) ModOde2(t,y,immunoreactivity, perTV, perNB,...
    kAR, k_AR, Vv, Vs, n, CL, S, D_T), unique(tspan),ini,options);


%*** ODEs ARE SOLVED. NOW, GET AUC[C_IAR], TR, etc ***

%Figure out correction for physical decay
kI = log(2)/t_half;
cI = cI0*exp(-kI * t); %cI and cIO are in units of mCi/mg


%Used to calculate AUC[C_IA]
%Convert from concentration to radiation, while applying physical decay
%correction (with cI). y is in mol/L, MM in Da, and cI in mCi/mg (Ci/g)
cV = y(:,1) + y(:,2); %Ventricle free antibody (with and without immunoreactivity), mol/L
cS = y(:,3) + y(:,4); %Subarachnoid free antibody (with and without immunoreactivity), mol/L
cIV = 1000*MM*cV.*cI; %ventricular "bad" radiation, mCi/L
cIS = 1000*MM*cS.*cI; %subarachnoid "bad" radiation, mCi/L
cIA = cIV + cIS; %total radiation
cIobs = cIV;


%Used to calculate AUC[C_IAR] or to return
%Convert from concentration to radiation, while applying physical decay
%correction (with cI). y is in mol/L, MM in Da, and cI in mCi/mg (Ci/g)
cIBV = 1000*MM*y(:,1).*cI;    %Free antibody w/  IR in ventricles; mCi/L
cINV = 1000*MM*y(:,2).*cI;    %Free antibody w/o IR in ventricles; mCi/L
cIBS = 1000*MM*y(:,3).*cI;    %Free antibody w/  IR in subarachnoid; mCi/L
cINS = 1000*MM*y(:,4).*cI;    %Free antibody w/o IR in subarachnoid; mCi/L
cIARV = 1000*MM*y(:,5).*cI;   % mCi/L
%RV doesn't need to be calculated
cIARS = 1000*MM*y(:,7).*cI;   % mCi/L
%RS doesn't need to be calculated
cINARV = 1000*MM*y(:,9).*cI;  % mCi/L
%NRV doesn't need to be calculated
cINARS = 1000*MM*y(:,11).*cI; % mCi/L
%NRS doesn't need to be calculated
returns = [cIBV cINV cIBS cINS cIARV cIARS cINARV cINARS]; % mCi/L

%Add the extra columns if needed. Saves about 3ms of time.
if return_all
    cIRV =  1000*MM*y(:,6).*cI;
    cIRS =  1000*MM*y(:,8).*cI;
    cINRV= 1000*MM*y(:,10).*cI;
    cINRS= 1000*MM*y(:,12).*cI;
    returns = [returns cIRV cIRS cINRV cINRS]; % mCi/L
end



%Find radiation values that are volume weighted, so just in units mCi
rIV = (cIBV+cINV)*Vv;
rIS = (cIBS + cINS)*Vs;
rIARV = (cIARV)*Vv;
rIARS = (cIARS)*Vs;
rIAR = rIARV+rIARS;
rIA = rIV+rIS;
AUC_rIAR = trapz(t,rIAR); %rIAR is in mCi and AUC is in mCi seconds
AUC_rIA = trapz(t,rIA);
rTR = AUC_rIAR/AUC_rIA;



if size(time_values_needed) ~= size(cIobs)
    cIobs = cIobs';
end


AUC_cIARV = trapz(t,cIARV); %cIARV is in mCi/L so AUC is in mCi s/L
AUC_cIARS = trapz(t,cIARS);
AUC_cIV = trapz(t,cIV);
AUC_cIS = trapz(t,cIS);

AUC_cIAR = (AUC_cIARV + AUC_cIARS);
AUC_cIA = (AUC_cIV + AUC_cIS);
TherapRatio =  AUC_cIAR/ AUC_cIA;


%Create an event so function stops when it's over, instead at a set time
function [pos,isterm,dir] = eventFcn(t,y, cDose) %cDose in mol/L
concentrations = y([1:5 7 9 11]);
%Pos is the value that we want to be zero
%0.001 is right to 5 digits compared to 0.00001
pos = max(0, sum(concentrations)-.001*cDose); %Ignore receptors at 6,8,10,12
isterm = 1;  % Is terminal, halt integration 
dir = 0;   % The zero can be approached from either direction

