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
% This file is similar to solveModel.m, except that it can account for
% multiple infusion times (and with varying activities administered). It
% essentially works by doing the actions of solveModel.m, but by iterating
% through all the infusion times and then adding up the outputs. This may
% not be mathematically exact (I haven't gone through the math to be sure,
% but intuitively, I'm not sure that the solution to the sum of ODEs is the
% same as the sum of the solutions to the ODEs; especially because each ODE
% assumes all receptors are unbound at the start), but it should yield
% approximately the same results because by the next infusion time (at
% least among the range of times I explored and are most clinically
% reasonable), most of the antigen receptors are already unbound.

% IMPORTANT- If you change this file, don't forget to change
% solveModel.m as well

% ARGUMENTS
% timeAdmins
%               This argument is the only one special to the
%               solveMultiInfus.m file (and not in the solveModel.m file).
%               This is a vector of time points (units: seconds) at which
%               the activity is administered. This is not in the
%               solveModel.m file because there activity is assumed to only
%               be administered at time t=0. This vector should be the same
%               length as the activityAdmins vector input.
% activityAdmins
%               Vector of amount of activity administered (units: mCi
%               [1mCi=37MBq]). Unlike solveModel.m, this is a vector and
%               the variable name is plural. Each element of the vector
%               should correlate with the respective indexed element in the
%               timeAdmins vector.
% Vv:           Ventricular volume (units: liters)
% CL            Clearance (L/s)
% time_values_needed
%               Vector of time values to run the model for and get the
%               associated outputs for to plot (units: seconds)
%               Note- This can actually change the AUC reported.
% defined_params
%               Vector of all other parameters. In order:
%   1: immunoreactivity
%               Immunoreactivity (unitless, on [0,1))
%   2: tumorload
%               Initial antigen concentration. Aka R0, but units: mol/L
%   3: perTV    TumorVent, percent of tumor in ventricles (unitless, on [0,1])
%   4: perNB    Percent nonspecificity (unitless, [0,1])
%   5: kAR      Association rate constant of forming antigen-receptor complex
%               (units: 1/(M s))
%   6: k_AR     Dissociation rate constant of forming antigen-receptor complex
%               (units: 1/s)
%   7: V        total CSF volume (units: liters)
%   8: n
%   9: MM       Molar mass (units: daltons)
%   10: cI0     Specific activity (units: mCi/mg, Ci/g)
%   11: t_half  Half-life (units: sec)
%   12: S       Optional. Surface area of entire CSF space (units: m^2)
%           
% X return_all  "return_all" isn't available in this file (didn't get
%               around to adding it) 



% OUTPUTS
% t:        Times outputted by the ODE solver. Should match
%           time_values_needed, but I haven't tested it thoroughly. (Units:
%           sec)
% cIobs:
%           Total radiation we can observe, i.e. that is only in the
%           ventricles. Equivalent to cIV. (unit: mCi/L)
% **SEE NOTATION DESCRIPTION BELOW TO UNDERSTAND REMAINING OUTPUTS**
%       AUC_f means the AUC[f] where AUC stands for the area under the
%       curve of function f. Lowercase c and r stand for concentration
%       (mol/L) or radiation (mCi), respectively. The cIx indicates that we
%       are talking about the radioactivity concentration of species X
%       (mCi/L). Thus AUC_cIAR is the AUC of the specific activity of
%       species AR (mCi s/L). A = antigen. R = receptor. Suffixes like V
%       and S added to species names can be used to represent only those
%       parts coming from the ventricles or subarachnoid space
%       respectively.
% AUC_cIAR:
% AUC_cIA:
% returns:
%           Many of the variables are not provided directly, but in a
%           "returns" vector variable because these values are not often
%           necessary. They are only provided for special analyses. The
%           contents are [cIBV cINV cIBS cINS cIARV cIARS cINARV cINARS]. If
%           returns_all is true, then cIRV, cIRS, cINRV, cINRS are also
%           included at the end. (units: mCi/mL)



% DETAILS
% Tested in MATLAB® 2017a for Mac.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [t, cIobs, AUC_cIAR, AUC_cIA, returns] = solveModelMultInfus(...
    timeAdmins, activityAdmins, ...
    Vv, CL, all_time_values_needed, ...
    defined_params) 

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
    tumorload = 4.755e-7; %aka R0
    perTV = 0; %TumorVent, percent of tumor in ventricles
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

% && has to be used, not & in order to use short-circuiting
if exist('defined_params', 'var') && length(defined_params) > 11
    S = defined_params(12) * 100; %convert from m2 to dm2
else
    %Surface Area of entire CSF space
    %Scaled knowing cited value that 140mL CSF volume has 0.18m2 (1800cm2)
    %Multiply by 100 to convert below units from m2 to dm2 so that units
    %become L when you multiply with D_T
    S = 100*0.18*sqrt((Vv+Vs)./0.14);
end

% Units in dm so becomes L when multiply with surface area (in dm2)
D_T = 0.0001; %Tumor cell spherical diameter (0.0001dm = 10um)


% Initial Values
%dose = (9.52e-8/2*0.140/Vv) * 3.84; %mol/L
cDose=activityAdmins(1)/cI0*1/1000*1/MM*1/Vv; %activityAdmin is in mCi
y = [0, 0, 0, 0, ...
    0, perTV * tumorload, 0, (1-perTV) * tumorload, ...
    0, perNB * (Vv/V) * tumorload, 0, perNB * (Vs/V) * tumorload];
t=0;
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

timeAdmins = [timeAdmins Inf]; % Add inf so the > and < time condition to isolate the time range still works
% Iterate through each of the administration times.
for i = 1:(length(timeAdmins)-1)
    
    % Inputs
    cDose=activityAdmins(i)/cI0*1/1000*1/MM*1/Vv ; %activityAdmin is in mCi
    ini = y(end, :);
    ini(1:2) = ini(1:2) + [immunoreactivity * cDose, (1-immunoreactivity) * cDose];
    %ini = [cBV0, cNV0, cBS0, cNS0, cARV0, cRV0, cARS0, cRS0, cNARV0, cNRV0, cNARS0, cNRS0];
    
    time_values_needed = all_time_values_needed(...
        all_time_values_needed(:) >= timeAdmins(i)...
        & all_time_values_needed(:) < timeAdmins(i+1));
    if time_values_needed(1) == timeAdmins(i)
        tspan = time_values_needed;
    elseif isnan(time_values_needed(1))
        %NaN start indicates ode isn't supposed to start from time 0 and thus
        %should bypass the check above
        tspan = [time_values_needed(2:end)];
    else
        tspan = [timeAdmins(i) time_values_needed];
    end
    if i < length(timeAdmins)-1 %less than so inf isn't added
        tspan = [tspan timeAdmins(i+1)];
    end
    %tolerance = 1e-16;
    tolerance = 2.2205e-14;
    
    %[t,y] = ode45(@(t,y) ModOde2(t,y,immunoreactivity, perTV, perNB,...
    %    kAR, k_AR, Vv, Vs, n, CL), unique(tspan),ini,options);
    
    if i==length(timeAdmins)-1 %if this is the last iteration
        %event stops the ode when values are low enough
        options = odeset('RelTol', tolerance, 'AbsTol', tolerance);
        %,...'Events', @(t,y)eventFcn(t,y,cDose)
    else
        options = odeset('RelTol', tolerance, 'AbsTol', tolerance);
    end
    
    myode = @(tm,yv) ModOde2(tm, yv,  immunoreactivity, perTV, perNB,...
        kAR, k_AR, Vv, Vs, n, CL, S, D_T);
        
    sol = ode45(myode,[min(tspan) max(tspan)],ini,options);
    tnew = tspan';
    ynew=deval(sol,tspan)';
    
    %[tnew,ynew] = ode113(@AbFinal,unique(tspan),ini,options, immunoreactivity, perTV, perNB, kAR, k_AR, Vv, Vs, n, CL);
    
    y = [y(1:end-1, :); ynew(1:end, :)];
    t = [t(1:end-1, :); tnew(1:end, :)];
end

% Figure out correction for physical decay
kI = log(2)/t_half;
cI = cI0*exp(-kI * t);

cV = y(:,1) + y(:,2); %mol/L
cS = y(:,3) + y(:,4); %mol/L

% Used to calculate AUC[C_IAR] or to return
% Convert from concentration to radiation, while applying physical decay
% correction (with cI). y is in mol/L, MM in Da, and cI in mCi/mg (Ci/g)
cIBV = 1000*MM*y(:,1).*cI;    %Free antibody w/  IR in ventricles; mCi/L
cINV = 1000*MM*y(:,2).*cI;    %Free antibody w/o IR in ventricles; mCi/L
cIBS = 1000*MM*y(:,3).*cI;    %Free antibody w/  IR in subarachnoid; mCi/L
cINS = 1000*MM*y(:,4).*cI;    %Free antibody w/o IR in subarachnoid; mCi/L
cIARV = 1000*MM*y(:,5).*cI;   % mCi/L
% RV doesn't need to be calculated
cIARS = 1000*MM*y(:,7).*cI;   % mCi/L
% RS doesn't need to be calculated
cINARV = 1000*MM*y(:,9).*cI;  % mCi/L
% NRV doesn't need to be calculated
cINARS = 1000*MM*y(:,11).*cI; % mCi/L
% NRS doesn't need to be calculated
returns = [cIBV cINV cIBS cINS cIARV cIARS cINARV cINARS]; % mCi/L

% Do something similar with combined values
cIV = 1000*MM*cV.*cI; %ventricular radiation, radiation observed
cIS = 1000*MM*cS.*cI;
cIVS = cIV + cIS; %total radiation
cIobs = cIV;

if size(time_values_needed) ~= size(cIobs)
    cIobs = cIobs';
end

%figure;
%plot(t/3600, cIV);
%hold on;
%plot(t/3600, cIS);
%legend(['cIV ' num2str(trapz(cIV))],['cIS ' num2str(trapz(t,cIS))]);

AUC_cIARV = trapz(t,cIARV);
AUC_cIARS = trapz(t,cIARS);
AUC_cIV = trapz(t,cIV);
AUC_cIS = trapz(t,cIS);

AUC_cIAR = (AUC_cIARV + AUC_cIARS);
AUC_cIA = (AUC_cIV + AUC_cIS);
TherapRatio =  AUC_cIAR/ AUC_cIA; %unitless


% Create an event so function stops when it's over, instead at a set time
function [pos,isterm,dir] = eventFcn(t,y, cDose) %cDose in mol/L
concentrations = y([1:5 7 9 11]);
%Pos is the value that we want to be zero
%0.001 is right to 5 digits compared to 0.00001
pos = max(0, sum(concentrations)-.001*cDose); %Ignore receptors at 6,8,10,12
isterm = 1;  % Is terminal, halt integration 
dir = 0;   % The zero can be approached from either direction
