function [validStage, powerExtracted, maxMach, degOfRxnHub, bladeCountStator, bladeCountRotor] = analyzeStageFunction(kv2,alpha2, lambda, phi)
%ANALYZESTAGE determines the performance metrics of a turbine stator-rotor
%stage.  Used to iterate.  see analyzeStage for Nomenclature, Flow
%convention and Unit Convention.  Stage Requirements are constant, defined
%in 'Stage Requirements' 
validStage = true;
%% Code Performance and Plotting Parameters
printValues = false;
printTables = false;
radRes = 3; %radial point count for equilibrium. recommend odd, >2
%axialRes = 3; %axial resolution for plots, points per blade

%% **** Design Constants to modify the performance of the stage ********
%kv1, alpha2, lambda and phi from input parameters
kr2 = .965; % determines r2 in % of r3
kr1 = 1; % determines r1 in % of r2

alpha1 = 0; %design constraint, first stage is axial flow

statorAR = .8; %.5 for hpt and 4-5 for lpt
rotorAR = 1.05;

statorTaper = 1.2;
rotorTaper = .87;

geoSpacing = 1.2*0.0254; %.75 inches recommended, converted to SI
zweifel = .8; % recommended starting point

%% Inlet Conditions [SI]
t01 = 1700;
p01 = 2050* 10^3;
M1 = .4;

%% Stage Requirements
eta_tt = .91;
stageTotalPressRatio = 3.2;
statorLossCoeff = .053;
rpm = 13000;
omega = (rpm*2*pi/60);
m_dot = 70;

%% Average Values (can iterate later if needed)
R = 260;
gamma = 1.33;
Cp = (gamma*R)/(gamma-1);

%% 1 - Isen to find t1,p1
p1 = p01/(1+(gamma-1)/2*M1^2)^(gamma/(gamma-1));
t1 = t01/(1+(gamma-1)/2*M1^2);

%% 2 - stator (1->2) Adiabatic Exp, no work
v1 = [M1*sqrt(gamma*R*t1),0]; %all axial;
A1 = m_dot*R*t1/(p1*v1(1));

%% 3 - Find absolute conditions @2 via Stator total pressure loss coeff
p02 = p01 - (statorLossCoeff*.5*p1/(R*t1) * norm(v1)^2);
t02 = t01;

%% 4/5 - Rotor 1 Thermodynamics (2->3) Adiabatic Exp, work extraction
p03 = p01/stageTotalPressRatio;

Tau0 = 1-eta_tt*(1-(p03/p01)^((gamma-1)/gamma));
t03 = Tau0*t01;

entropyChange = (t01-t03)*Cp;

%% 7/8/9a - Um3 via lambda,  v3 via phi and r3 via u3
u3 = [0,-sqrt(entropyChange/lambda)]; 
v3 = [norm(u3)*phi, -inf];  % TO-DO: this involves phi - can this be negative or positive?

v2 = [(v3(1)-v1(1))*kv2 + v1(1), -inf]; %kv2 is design parameter
v2(2) = v2(1)*tan(alpha2);

%% 9b - rotor constants and thus rm2,rm1
rm3 = norm(u3)/omega;

rm2 = rm3 * kr2;
rm1 = rm2 * kr1;

%% 10 - um2 from omega and the radius
u2 = [0,-omega*rm2];
u1 = [0,0] ;% stationary

%% 12/14 - finish velocity triangles for the rotor

v3(2) = (entropyChange+(-u2(2)*v2(2)))/-u3(2); %Euler Turbo Eqn, neg u due to conventions
entCheck = u2(2)*v2(2) - u3(2)*v3(2);

if(abs(entCheck -entropyChange) > 1) % leave small room for rounding/truncating
    disp("ERROR: Due to unconventional angles/geometry Euler Turbo Eqn has produced an invalid results as implemented. Results invalid");
end

w2 = v2-u2;
w3 = v3-u3;

beta2 = atan(w2(2)/w2(1));
beta3 = atan(w3(2)/w3(1));
alpha3 = atan(v3(2)/v3(1));

%% 11/15 - determine static conditions into and out of rotor (2) & (3) 

t2 = t02 - .5* norm(v2)^2/Cp; % TO-DO: could these be the problem?
M2 = norm(v2)/sqrt(gamma*R*t2);
p2 = p02 / (1+(gamma-1)/2*M2^2)^(gamma/(gamma-1));
A2 = m_dot*R*t2/(norm(v2)*p2);

t3 = t03 - .5* norm(v3)^2/Cp;  % TO-DO: this is off by about 5% due to Cp not varying with termperature (or gamma)
M3 = norm(v3)/sqrt(gamma*R*t3);
p3 = p03 / (1+(gamma-1)/2*M3^2)^(gamma/(gamma-1));
A3 = m_dot*R*t3/(norm(v3)*p3);%% 7/8/9a - Um3 via lambda,  v3 via phi and r3 via u3
u3 = [0,-sqrt(enthalpyChange/lambda)]; 
v3 = [norm(u3)*phi, -inf];  % TO-DO: this involves phi - can this be negative or positive?

v2 = [(v3(1)-v1(1))*kv2 + v1(1), -inf]; %kv2 is design parameter
v2(2) = v2(1)*tan(alpha2);

%% 9b - rotor constants and thus rm2,rm1
rm3 = norm(u3)/omega;

rm2 = rm3 * kr2;
rm1 = rm2 * kr1;

%% 10 - um2 from omega and the radius
u2 = [0,-omega*rm2];
u1 = [0,0] ;% stationary

%% 12/14 - finish velocity triangles for the rotor TO-DO: This still seems incorrect with dot products

v3(2) = (enthalpyChange -(u2(2)*v2(2)))/-u3(2); %Euler Turbo Eqn, neg u due to conventions
entCheck = u2(2)*v2(2) -u3(2)*v3(2);

if(entCheck ~=enthalpyChange)
    disp("ERROR: Due to unconventional angles/geometry Euler Turbo Eqn has produced an invalid results as implemented. Results invalid");
end

w2 = v2-u2;
w3 = v3-u3;

beta2 = atan(w2(2)/w2(1));
beta3 = atan(w3(2)/w3(1));
alpha3 = atan(v3(2)/v3(1));

%% 11/15 - determine static conditions into and out of rotor (2) & (3) 

t2 = t02 - .5* norm(v2)^2/Cp; % TO-DO: could these be the problem?
M2 = norm(v2)/sqrt(gamma*R*t2);
p2 = p02 / (1+(gamma-1)/2*M2^2)^(gamma/(gamma-1));
A2 = m_dot*R*t2/(v2(1)*p2);

t3 = t03 - .5* norm(v3)^2/Cp;  % TO-DO: this is off by about 5% due to Cp not varying with termperature (or gamma)
M3 = norm(v3)/sqrt(gamma*R*t3);
p3 = p03 / (1+(gamma-1)/2*M3^2)^(gamma/(gamma-1));
A3 = m_dot*R*t3/(v3(1)*p3);

%% 16 - Geometry - find the width of the channels given areas and mid radius
rh1 = rm1 - A1/(4*pi*rm1);
rt1 = 2*rm1 - rh1;

rh2 = rm2 - A2/(4*pi*rm2);
rt2 = 2*rm2 - rh2;

rh3 = rm3 - A3/(4*pi*rm3);
rt3 = 2*rm3 - rh3;

%% 16b - find velocity triangles as function of radius
r1_fun = linspace(rh1,rt2,radRes);
rStator_fun = linspace(mean([rh1,rh2]),mean([rt1,rt2]),radRes);
r2_fun = linspace(rh2,rt2,radRes);
rRotor_fun = linspace(mean([rh2,rh3]),mean([rt2,rt3]),radRes);
r3_fun = linspace(rh3,rt3,radRes);

% v<_>_fun adn u<_>fun via Free-Vortex
v1_fun = [v1(1)*ones(1,radRes); rm1*v1(2)./r1_fun];
v2_fun = [v2(1)*ones(1,radRes); rm2*v2(2)./r2_fun];
v3_fun = [v3(1)*ones(1,radRes); rm3*v3(2)./r3_fun];

u2_fun = [zeros(1,radRes);  -omega*r2_fun];
u3_fun = [zeros(1,radRes); -omega*r3_fun];

w2_fun = v2_fun - u2_fun;
w3_fun = v3_fun - u3_fun;

%beta1_fun = zeros(size(w2_fun(1,:)));
beta2_fun = atan(w2_fun(2,:)./w2_fun(1,:));
beta3_fun = atan(w3_fun(2,:)./w3_fun(1,:));
alpha1_fun = atan(v1_fun(2,:)./v1_fun(1,:)); % design contraint - should be 0
alpha2_fun = atan(v2_fun(2,:)./v2_fun(1,:)); % <-alpha2 was set for design @ midradius
%alpha3_fun = atan(v3_fun(2,:)./v3_fun(1,:));

%% 17 - degree of reaction - mid and function - TO-DO: Ask attia if it is expected that the reacion at the tip is smaller than hub
% degOfRxnNumerator = .5*norm(w3)^2 - .5*norm(w2)^2;
% degOfRxnDenominator = .5*norm(v2)^2 - .5*norm(v3)^2 + .5*norm(w3)^2 - .5*norm(w2)^2;
% degOfRxn = degOfRxnNumerator/degOfRxnDenominator;

degOfRxnNumerator_fun = .5*vecnorm(w3_fun).^2 - .5*vecnorm(w2_fun).^2;
degOfRxnDenominator_fun = .5*vecnorm(v2_fun).^2 - .5*vecnorm(v3_fun).^2 + .5*vecnorm(w3_fun).^2 - .5*vecnorm(w2_fun).^2;
degOfRxn_fun = degOfRxnNumerator_fun./degOfRxnDenominator_fun;

degOfRxnHub = degOfRxn_fun(1);

%% 18a - meridonal view of stage
bStator = rStator_fun(end)-rStator_fun(1);
bRotor = rRotor_fun(end)-rRotor_fun(1);

staggerStator_fun = (alpha1_fun+alpha2_fun)/2;
staggerRotor_fun = (beta2_fun+beta3_fun)/2;

% for reference - TO-DO: the chord seems too long for the given parameters
% cStatorMean = statorAR*bStator;
% cRotorMean = rotorAR*bRotor;

% root and tip are physical in this section, violate nomenclature for clarity
cAxialStatorRoot = (bStator/statorAR)./((1+(statorTaper-1)/2))*cos(staggerStator_fun(end)); 
cAxialStatorTip = cAxialStatorRoot*statorTaper;
cAxialStator_fun = linspace(cAxialStatorTip,cAxialStatorRoot,radRes);
%cStator_fun = cAxialStator_fun./cos(staggerStator_fun);

cAxialRotorRoot = (bRotor/rotorAR)./(1+(rotorTaper-1)/2)*cos(staggerRotor_fun(1));
cAxialRotorTip = cAxialRotorRoot*rotorTaper;
cAxialRotor_fun = linspace(cAxialRotorRoot,cAxialRotorTip,radRes);
%cRotor_fun = cAxialRotor_fun./cos(staggerRotor_fun);

% statorTaper = cAxialStator_fun(1)/cAxialStator_fun(end); %inv b/c radius convention. tip/root
% rotorTaper = cAxialRotor_fun(end)/cAxialRotor_fun(1);

%% 18d - plot the position of root, half and tip with 2nd order polynomial - TO-DO: Finish this with LE&TE points
%radiusPloty  = [r1_fun; rStator_fun; rStator_fun;r2_fun; rRotor_fun;rRotor_fun;r3_fun];
radiusPlotxIncruments = [0, geoSpacing/2-mean(cAxialStator_fun)/2, mean(cAxialStator_fun), geoSpacing/2-mean(cAxialStator_fun)/2 ,geoSpacing/2 - mean(cAxialRotor_fun)/2,mean(cAxialRotor_fun), geoSpacing/2 - mean(cAxialRotor_fun)/2];
radiusPlotx = cumsum(radiusPlotxIncruments);

radiusToInterpY = [r1_fun; rStator_fun; r2_fun; rRotor_fun; r3_fun];
radiusToInterpXInc = [0,geoSpacing/2,geoSpacing/2,geoSpacing/2,geoSpacing/2];
radiusToInterpX = cumsum(radiusToInterpXInc);

%xRange = linspace(0,max(radiusPlotx),100);
newRadii = zeros(length(radiusPlotx),radRes);

for i=1:radRes
    polyRad = polyfit(radiusToInterpX,rot90(radiusToInterpY(:,i)),2); %2nd order poly fit
    newRadii(:,i) = polyval(polyRad, radiusPlotx);
    %yVals = polyval(polyRad, xRange);
end

% Interpolated Meridional View
statorYInterp = [rot90(newRadii(2,:),2),rot90(fliplr(newRadii(3,:)),2)];
rotorYInterp = [newRadii(5,:),fliplr(newRadii(6,:))];

%% 18e - finalize airfoils per row
sRotor = zweifel * mean(cAxialRotor_fun) / (2*cos(beta3)^2 * (tan(beta2)-tan(beta3))); % TO-DO: determine convention that zweifel used for trig (pos or neg angles)
sStator = zweifel * mean(cAxialStator_fun) / (2*cos(alpha2)^2 * (tan(alpha1)-tan(alpha2)));

indexOfLEMidSpan = ceil(radRes/2);
rRotorLEMidSpan = rotorYInterp(indexOfLEMidSpan);
rStatorLEMidSpan = statorYInterp(indexOfLEMidSpan);

bladeCountRotor = floor(abs(2*pi*rRotorLEMidSpan/sRotor));
bladeCountStator = floor(abs(2*pi*rStatorLEMidSpan/sStator));

%% Misc calculations for design constraint validation
Mw2 = norm(w2)/sqrt(gamma*R*t2); %relative mach numbers
Mw3 = norm(w3)/sqrt(gamma*R*t3);

powerExtracted = dot(u3,v3)+dot(u2,v2);

%% Check for design constraints
if(abs(degOfRxn_fun(1)-.1)>.01)
    validStage = false;
    %fprintf("POTENTIAL DESIGN FLAW! your Degree of reaction at the hub is %.3f%%, which is outside of the required range of 10%% +- 1%%\n",degOfRxn_fun(1)*100);
end

if(abs(statorTaper-1.2)>.1)
       % validStage = false;
   %fprintf("POTENTIAL DESIGN FLAW! your stator taper ratio is %.4f, which is outside of the required range of 1.1-1.3\n",statorTaper); 
end

if(abs(rotorTaper-.85)>.05)
        %validStage = false;
    %fprintf("POTENTIAL DESIGN FLAW! your rotor taper ratio is %.4f, which is outside of the required range of .8-.9\n",rotorTaper);
end

if(abs(statorAR-1.05)>.25)
        %validStage = false;
    %fprintf("POTENTIAL DESIGN FLAW! your stator aspect ratio is %.3f, which is outside of the required range of .8-1.3\n",statorAR);
end

if(abs(rotorAR-1.05)>.25)
        %validStage = false;
    %fprintf("POTENTIAL DESIGN FLAW! your rotor aspect ratio is %.3f, which is outside of the required range of .8-1.3\n",rotorAR);
end

if(abs(phi-.95)>.151)
        validStage = false;
    fprintf("POTENTIAL DESIGN FLAW! your flow coeficcient (phi) is %.2f, which is outside of the required range of .8-1.1\n",phi);
end

%% Check for recommended constraints (from the notes)
maxMach = max([M1,M2,M3,Mw2,Mw3]);
if(maxMach > 1.34) %TBD
        validStage = false;
   %fprintf("POTENTIAL DESIGN FLAW! your maximum local or abs Mach is %.2f, which is above the suggested M = 1.3 cap of the notes\n",maxMach);
end

maxAngle = max([beta2,beta3,alpha1,alpha2,alpha3]);
if(rad2deg(maxAngle)>72)
        validStage = false;
    %fprintf("POTENTIAL DESIGN FLAW! The highest angle is %.2f [deg], which is above the suggested 72 [deg] of the notes\n",rad2deg(maxAngle));
end

statorTurnAng = abs(alpha2-alpha1);
rotorTurnAng = abs(beta3-beta2);
if(statorTurnAng>deg2rad(120))
        validStage = false;
    %fprintf("POTENTIAL DESIGN FLAW! your stator turn angle is %.1f, which is above the suggested 120 [deg] max of the notes\n",rad2deg(statorTurnAng));
end

if(rotorTurnAng>deg2rad(120))
        validStage = false;
    %fprintf("POTENTIAL DESIGN FLAW! your rotor turn angle is %.1f, which is above the suggested 120 [deg] max of the notes\n",rad2deg(rotorTurnAng));
end

if(abs(bladeCountRotor-bladeCountStator)<4)%TBD: CAN MODIFY THIS WITH  RADII LATER< LETTING IT BE NOW
       % validStage = false;
    %fprintf("POTENTIAL DESIGN FLAW! you have  %.0f stator blades and %.0f rotor blades.  Recommend having at least 4 blades difference for harmonics\n",bladeCountStator, bladeCountRotor);
end

%% Check for extreme errors
if(~isreal(A2) || ~isreal(A3))
    validStage = false;
end

minTemp = min([t1,t2,t3]);
if (minTemp < 0)
    validStage = false;
end

%% Print Values if printValues
if printValues
    fprintf("Total Pressure Values [Pa]:\n\t\tp01\t\t\t\tp02\t\t\t\tp03\n\t\t%.4e\t\t%.4e\t\t%.4e\n",p01,p02,p03);
    fprintf("Static Pressure Values [Pa]:\n\t\tp1\t\t\t\tp2\t\t\t\tp3\n\t\t%.4e\t\t%.4e\t\t%.4e\n",p1,p2,p3);
    fprintf("Total Temp Values [K]:\n\t\tt01\t\t\t\tt02\t\t\t\tt03\n\t\t%.1f\t\t\t%.1f\t\t\t%.1f\n",t01,t02,t03);
    fprintf("Static Temp Values [K]:\n\t\tt1\t\t\t\tt2\t\t\t\tt3\n\t\t%.1f\t\t\t%.1f\t\t\t%.1f\n",t1,t2,t3);    
    fprintf("Deflection Angles,alpha [rad] :\n\t\talpha_1\t\t\talpha_2\t\t\talpha_3\n\t\t%.3f\t\t\t%.3f\t\t\t%.3f\n",alpha1,alpha2,alpha3);
    fprintf("Relative Angles,beta [rad]:\n\t\tbeta_1\t\t\tbeta_2\t\t\tbeta_3\n\t\t-N/A-\t\t\t%.3f\t\t\t%.3f\n",beta2,beta2);
    fprintf("Degree of Reaction:\n\t\t%.2f%%\n",degOfRxn*100);
    fprintf("Flow Coeff:\n\t\t%.2f\n",phi);
    fprintf("Work Coeff:\n\t\t%.2f\n",lambda);
    fprintf("Mis Span Zweifel Coeff:\n\t\t%.2f\n",zweifel);
    fprintf("Stator Blade Count:\n\t\t%.0f\n",bladeCountStator);
    fprintf("Rotor Blade Count:\n\t\t%.0f\n",bladeCountRotor);
end

%% Print Tables (colums so they can be copy/pasted into spreadsheet/word
if printTables
    
fprintf("First Table with design choices:\n");
fprintf("\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n",rpm, zweifel, phi,v2(1)/v1(1)*100,v3(1)/v2(1)*100, rm3/rm2, rm2/rm1, rad2deg(alpha3), rotorAR, rotorTaper, statorAR, statorTaper);

%fprintf("Third Table with Velocity Triangle Info\n");
%fprintf("\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f",rad2deg(alpha1),-inf,u1(2),v1(1),v1(2),norm(v1),w1(1),w1(2),norm(w1));

end

%% Return Important Variables - Don't think this is required here
return;
end

