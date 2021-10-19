% Matthew Liepke, AE 403 Final Design Project initial code to validate some
% relations and to hopefully get a good understanding of what is required
% for this project

% FLOW CONVENTION: [axial, circumferential] or [Vax, Vu] for all velocities
% are avt/mid span.  The rotor spins in the neg(-) Vu direction, meaning
% the relative frame of the rotor moves at a neg(-) circumferential
% velocity

% UNIT CONVENTION: SI, absolute, radians.  This means K, Pa, m, radians, ect.

% NOMENCLATURE CONVENTION:
%   * (1) is pre-stator; (2) is post-stator, pre-rotor; (3) is post-rotor
%   * unless otherwise noted, values are average
%   * _fun are variables that vary over a given input (radius),
%       from low value (root) to high value (tip)
%   * tip radius (r_t) is greater than root radius (r_r) even though it may contradict a physical setup

clear;clc;close all;

%% Code Performance and Plotting Parameters
printValues = true;
printTables = true;
radRes = 5; %radial point count for equilibrium. recommend odd, >2
axialRes = 15; %axial resolution for plots, points per blade

%% **** Design Constants to modify the performance of the stage ********

kv2 = .5051; % determines v2ax in % between v1ax and v3ax
kr2 = 1.02034; % determines r2 in % of r3
kr1 = 1.02952; % determines r1 in % of r2

alpha1 = 0; %design constraint, first stage is axial flow
alpha2 = -1.18239; % guesswork value, angle out of stator and into rotor

statorAR = .8; %.5 for hpt and 4-5 for lpt
rotorAR = 1.05;

statorTaper = 1.2;
rotorTaper = .87;

geoSpacing = 2.8*0.0254; %.75 inches recommended, converted to SI
zweifel = .8; % recommended starting point

lambda = 2.06124; % work coeff, max is 2.2, can be lower if more eff needed
phi = .80604; %flow coeff, .8-1.1

%% Inlet Conditions [SI]
t01 = 1700;
p01 = 2050* 10^3;
M1 = .4;

%% Stage Requirements
eta_tt = .91;
stageTotalPressRatio = 3.2;
statorLossCoeff = .053;
rpm = 12750;
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

enthalpyChange = (t01-t03)*Cp;

%% 7/8/9a - Um3 via lambda,  v3 via phi and r3 via u3
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

t3 = t03 - .5* norm(v3)^2/Cp;  
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
r1_fun = linspace(rh1,rt1,radRes);
rStator_fun = linspace(mean([rh1,rh2]),mean([rt1,rt2]),radRes);
r2_fun = linspace(rh2,rt2,radRes);
rRotor_fun = linspace(mean([rh2,rh3]),mean([rt2,rt3]),radRes);
r3_fun = linspace(rh3,rt3,radRes);

% v<_>_fun adn u<_>fun via Free-Vortex
v1_fun = [v1(1)*ones(1,radRes); rm1*v1(2)./r1_fun];
v2_fun = [v2(1)*ones(1,radRes); rm2*v2(2)./r2_fun];
v3_fun = [v3(1)*ones(1,radRes); rm3*v3(2)./r3_fun];

u1_fun = [zeros(1,radRes); zeros(1,radRes)];
u2_fun = [zeros(1,radRes);  -omega*r2_fun];
u3_fun = [zeros(1,radRes); -omega*r3_fun];

w1_fun = [zeros(1,radRes); zeros(1,radRes)];
w2_fun = v2_fun - u2_fun;
w3_fun = v3_fun - u3_fun;

beta1_fun = zeros(size(w2_fun(1,:)));
beta2_fun = atan(w2_fun(2,:)./w2_fun(1,:));
beta3_fun = atan(w3_fun(2,:)./w3_fun(1,:));
alpha1_fun = atan(v1_fun(2,:)./v1_fun(1,:)); % design contraint - should be 0
alpha2_fun = atan(v2_fun(2,:)./v2_fun(1,:)); % <-alpha2 was set for design @ midradius
alpha3_fun = atan(v3_fun(2,:)./v3_fun(1,:));

%% 17 - degree of reaction - mid and function - TO-DO: DOR may not properly account for negative angles/velocities
degOfRxnNumerator = .5*norm(w3)^2 - .5*norm(w2)^2;
degOfRxnDenominator = .5*norm(v2)^2 - .5*norm(v3)^2 + .5*norm(w3)^2 - .5*norm(w2)^2;
degOfRxn = degOfRxnNumerator/degOfRxnDenominator;

degOfRxnNumerator_fun = .5*vecnorm(w3_fun).^2 - .5*vecnorm(w2_fun).^2;
degOfRxnDenominator_fun = .5*vecnorm(v2_fun).^2 - .5*vecnorm(v3_fun).^2 + .5*vecnorm(w3_fun).^2 - .5*vecnorm(w2_fun).^2;
degOfRxn_fun = degOfRxnNumerator_fun./degOfRxnDenominator_fun;

%% 18a - meridonal view of stage
bStator = rStator_fun(end)-rStator_fun(1);
bRotor = rRotor_fun(end)-rRotor_fun(1);

staggerStator_fun = (alpha1_fun+alpha2_fun)/2;
staggerRotor_fun = (beta2_fun+beta3_fun)/2;

% for reference - TO-DO: the chord seems too long for the given parameters
cStatorMid = statorAR*bStator;
cRotorMid = rotorAR*bRotor;

% root and tip are physical in this section, violate nomenclature for clarity
cAxialStatorRoot = (bStator/statorAR)./((1+(statorTaper-1)/2))*cos(staggerStator_fun(end)); 
cAxialStatorTip = cAxialStatorRoot*statorTaper;
cAxialStator_fun = linspace(cAxialStatorTip,cAxialStatorRoot,radRes);
cStator_fun = cAxialStator_fun./cos(staggerStator_fun);

cAxialRotorRoot = (bRotor/rotorAR)./(1+(rotorTaper-1)/2)*cos(staggerRotor_fun(1));
cAxialRotorTip = cAxialRotorRoot*rotorTaper;
cAxialRotor_fun = linspace(cAxialRotorRoot,cAxialRotorTip,radRes);
cRotor_fun = cAxialRotor_fun./cos(staggerRotor_fun);

statorTaper = cAxialStator_fun(1)/cAxialStator_fun(end); %inv b/c radius convention. tip/root
rotorTaper = cAxialRotor_fun(end)/cAxialRotor_fun(1);

%% 18d - plot the position of root, half and tip with 2nd order polynomial
radiusPloty  = [r1_fun; rStator_fun; rStator_fun;r2_fun; rRotor_fun;rRotor_fun;r3_fun];
radiusPlotxIncruments = [0, geoSpacing/2-mean(cAxialStator_fun)/2, mean(cAxialStator_fun), geoSpacing/2-mean(cAxialStator_fun)/2 ,geoSpacing/2 - mean(cAxialRotor_fun)/2,mean(cAxialRotor_fun), geoSpacing/2 - mean(cAxialRotor_fun)/2];
radiusPlotx = cumsum(radiusPlotxIncruments);

radiusToInterpY = [r1_fun; r2_fun; r3_fun];
radiusToInterpXInc = [0,geoSpacing,geoSpacing];
radiusToInterpX = cumsum(radiusToInterpXInc);

figure('Name','Interpolated Radius');
hold on;
axis equal;
plot(radiusToInterpX, radiusToInterpY);
xRange = linspace(0,max(radiusPlotx),100);

newRadii = zeros(length(radiusPlotx),radRes);
for i=1:radRes
    polyRad = polyfit(radiusToInterpX,rot90(radiusToInterpY(:,i)),2); %2nd order poly fit
    newRadii(:,i) = polyval(polyRad, radiusPlotx);
    yVals = polyval(polyRad, xRange);
    plot(radiusPlotx,newRadii(:,i),'o');
    plot(xRange, yVals, '--');
end

% UnInterpolated Meridional View
statorX = [geoSpacing/2 - cAxialStator_fun/2,fliplr(geoSpacing/2 + cAxialStator_fun/2)];
statorY = [rStator_fun, fliplr(rStator_fun)];

rotorX = [geoSpacing*3/2 - cAxialRotor_fun/2,fliplr(geoSpacing*3/2 + cAxialRotor_fun/2)];
rotorY = [ rRotor_fun, fliplr(rRotor_fun)];
    
figure('Name','Stage Meridional View');
hold on;
plot(radiusPlotx,radiusPloty);
plot(statorX, statorY, '--');
plot(rotorX, rotorY, '--');

title('Meridional View of Turbine Stage');
ylabel('radial distance [m]');
xlabel('axial distance [m]');
axis equal;
ylim([min(min(radiusPloty)),max(max(radiusPloty))]);
xlim([0,radiusPlotx(end)]); 

% Interpolated Meridional View
statorYInterp = [rot90(newRadii(2,:),2),rot90(fliplr(newRadii(3,:)),2)];
rotorYInterp = [newRadii(5,:),fliplr(newRadii(6,:))];

figure('Name','Interpolated Meridional View');
hold on;
plot(radiusPlotx,newRadii);
plot(statorX,statorYInterp,'--','Color','k');
plot(rotorX, rotorYInterp,'--','Color','k');

ylabel('radial distance [m]');
xlabel('axial distance [m]');
axis equal;
% ylim([0,.38]);
% xlim([0,.05]);

%% Determine Annulus Area at each points (including LE and TE of each blade)
A_fun = pi*(newRadii(:,end).^2 - newRadii(:,1).^2);

%% 18e - finalize airfoils per row
sRotor = zweifel * mean(cAxialRotor_fun) / (2*cos(beta3)^2 * (tan(beta2)-tan(beta3))); % TO-DO: determine convention that zweifel used for trig (pos or neg angles)
sStator = zweifel * mean(cAxialStator_fun) / (2*cos(alpha2)^2 * (tan(alpha1)-tan(alpha2)));

indexOfLEMidSpan = ceil(radRes/2);
rRotorLEMidSpan = rotorYInterp(indexOfLEMidSpan);
rStatorLEMidSpan = statorYInterp(indexOfLEMidSpan);

bladeCountRotor = ceil(abs(2*pi*rRotorLEMidSpan/sRotor))
bladeCountStator = ceil(abs(2*pi*rStatorLEMidSpan/sStator))

%% Extras - Plot Stator & Rotor Blade Shape (stagger with random airfoil)
% [x y z] = [axial (cAx), __ , radial]

staggerRotor_fun3D = zeros(radRes,axialRes);
z_rotor = zeros(radRes,axialRes);
x_rotor = zeros(radRes,axialRes);
y_rotor = zeros(radRes,axialRes);

staggerStator_fun3D = zeros(radRes,axialRes);
z_stator = zeros(radRes,axialRes);
x_stator = zeros(radRes,axialRes);
y_stator = zeros(radRes,axialRes);

for i=1:radRes
   
   staggerRotor_fun3D(i,:) = linspace(beta2_fun(i),beta3_fun(i),axialRes);
   zR = cumsum(tan(staggerRotor_fun3D(i,:) * cAxialRotor_fun(i)/axialRes));
   z_rotor(i,:) = zR;
   x_rotor(i,:) = linspace(-cAxialRotor_fun(i)/2,cAxialRotor_fun(i)/2,axialRes);
   y_rotor(i,:) = linspace(rotorYInterp(i),rotorYInterp(end-i+1),axialRes);

   staggerStator_fun3D(i,:) = linspace(beta1_fun(i),beta2_fun(i),axialRes);
   zS = cumsum(tan(staggerStator_fun3D(i,:) * cAxialStator_fun(i)/axialRes));
   z_stator(i,:) = zS;
   x_stator(i,:) = linspace(-cAxialStator_fun(i)/2,cAxialStator_fun(i)/2,axialRes);
   y_stator(i,:) = linspace(statorYInterp(i),statorYInterp(end-i+1),axialRes);


end

figure("Name","Rotor 3D Plot");
hold on;
surface(x_rotor,y_rotor,z_rotor,'edgecolor','k');

xlabel("Axial Direction");
ylabel("Radial Direction");
pbaspect([1 1 1]);

figure("Name","Stator 3D Plot");
hold on;
surface(x_stator,y_stator,z_stator,'edgecolor','k');

xlabel("Axial Direction");
ylabel("Radial Direction");
pbaspect([1 1 1]);

%% Misc calculations for design constraint validation
powerExtracted = m_dot*enthalpyChange;
fprintf("This stage currently extracts %2.f [W]\n",powerExtracted);

%% Rotor Relative Frame calculations
Mw2 = norm(w2)/sqrt(gamma*R*t2);
Mw3 = norm(w3)/sqrt(gamma*R*t3);

h02r = Cp*t2 + .5 * norm(w2)^2; 
h03r = Cp*t3 + .5 * norm(w3)^2;
t02r = h02r/Cp;
t03r = h03r/Cp;

p02r = p2 * (1 + (gamma-1)/2 * Mw2^2)^(gamma/(gamma-1));
p03r = p3 * (1 + (gamma-1)/2 * Mw3^2)^(gamma/(gamma-1));

%% Check for design constraints
if(abs(degOfRxn_fun(1)-.1)>.01)
    fprintf("POTENTIAL DESIGN FLAW! your Degree of reaction at the hub is %.2f%%, which is outside of the required range of 10%% +- 1%%\n",degOfRxn_fun(1)*100);
end

if(abs(statorTaper-1.2)>.1)
   fprintf("POTENTIAL DESIGN FLAW! your stator taper ratio is %.4f, which is outside of the required range of 1.1-1.3\n",statorTaper); 
end

if(abs(rotorTaper-.85)>.05)
    fprintf("POTENTIAL DESIGN FLAW! your rotor taper ratio is %.4f, which is outside of the required range of .8-.9\n",rotorTaper);
end

if(abs(statorAR-1.05)>.25)
    fprintf("POTENTIAL DESIGN FLAW! your stator aspect ratio is %.3f, which is outside of the required range of .8-1.3\n",statorAR);
end

if(abs(rotorAR-1.05)>.25)
    fprintf("POTENTIAL DESIGN FLAW! your rotor aspect ratio is %.3f, which is outside of the required range of .8-1.3\n",rotorAR);
end

if(abs(phi-.95)>.1501)
    fprintf("POTENTIAL DESIGN FLAW! your flow coeficcient (phi) is %.2f, which is outside of the required range of .8-1.1\n",phi);
end

%% Check for recommended constraints (from the notes)
maxMach = max([M1,M2,M3,Mw2,Mw3]);
if(maxMach > 1.3)
   fprintf("POTENTIAL DESIGN FLAW! your maximum local or abs Mach is %.2f, which is above the suggested M = 1.3 cap of the notes\n",maxMach);
end

maxAngle = max([beta2,beta3,alpha1,alpha2,alpha3]);
if(rad2deg(maxAngle)>72)
    fprintf("POTENTIAL DESIGN FLAW! The highest angle is %.2f [deg], which is above the suggested 72 [deg] of the notes\n",rad2deg(maxAngle));
end

statorTurnAng = abs(alpha2-alpha1);
rotorTurnAng = abs(beta3-beta2);
if(statorTurnAng>deg2rad(120))
    fprintf("POTENTIAL DESIGN FLAW! your stator turn angle is %.1f, which is above the suggested 120 [deg] max of the notes\n",rad2deg(statorTurnAng));
end
if(rotorTurnAng>deg2rad(120))
    fprintf("POTENTIAL DESIGN FLAW! your rotor turn angle is %.1f, which is above the suggested 120 [deg] max of the notes\n",rad2deg(rotorTurnAng));
end

if(abs(bladeCountRotor-bladeCountStator)<4)
    fprintf("POTENTIAL DESIGN FLAW! you have  %.0f stator blades and %.0f rotor blades.  Recommend having at least 4 blades difference for harmonics\n",bladeCountStator, bladeCountRotor);
end

%% Print Values if printValues - some of these only work correctly if radRes is 3
if printValues
    % H-S Stuff
    fprintf("Total Enthalpy [J / kg]:\n\t\th01\t\t\t\th02\t\t\t\th03\n\t\t%.4e\t\t%.4e\t\t%.4e\n",Cp*t01, Cp*t02,Cp*t03);
    fprintf("Static Enthalpy [J / kg]:\n\t\th1\t\t\t\th2\t\t\t\th3\n\t\t%.4e\t\t%.4e\t\t%.4e\n",Cp*t1, Cp*t2,Cp*t3);
    fprintf("Total Pressure Values [Pa]:\n\t\tp01\t\t\t\tp02\t\t\t\tp03\n\t\t%.4e\t\t%.4e\t\t%.4e\n",p01,p02,p03);
    fprintf("Static Pressure Values [Pa]:\n\t\tp1\t\t\t\tp2\t\t\t\tp3\n\t\t%.4e\t\t%.4e\t\t%.4e\n",p1,p2,p3);
    fprintf("Total Temp Values [K]:\n\t\tt01\t\t\t\tt02\t\t\t\tt03\n\t\t%.1f\t\t\t%.1f\t\t\t%.1f\n",t01,t02,t03);
    fprintf("Static Temp Values [K]:\n\t\tt1\t\t\t\tt2\t\t\t\tt3\n\t\t%.1f\t\t\t%.1f\t\t\t%.1f\n",t1,t2,t3);    
    fprintf("Entropy Change  Values, static [K]:\n\t\t 1-> 2 \t\t\t\t2 -> 3\n\t\t%.3f\t\t\t\t%.3f\n",Cp*log(t2/t1) - R*log(p2/p1),Cp*log(t3/t2) - R*log(p3/p2));
    fprintf("Entropy Change  Values, total [K]:\n\t\t 1-> 2 \t\t\t\t2 -> 3\n\t\t%.3f\t\t\t\t%.3f\n",Cp*log(t02/t01) - R*log(p02/p01),Cp*log(t03/t02) - R*log(p03/p02));
    %H-S, Relative Rotor Frame
    fprintf("Total Enthalpy, Relative Rotor [J / kg]:\n\t\th01r\t\t\t\th02r\t\t\t\th03r\n\t\t%.4e\t\t%.4e\t\t%.4e\n",nan,h02r, h03r);
    fprintf("Total Temp Values, Relative Rotor[K]:\n\t\tt01r\t\t\t\tt02r\t\t\t\tt03r\n\t\t%.4f\t\t\t%.4f\t\t\t%.1f\n",nan,t02r,t03r);
    fprintf("Total Pressure Values, Relative Rotor [Pa]:\n\t\tp01r\t\t\t\tp02r\t\t\t\tp03r\n\t\t%.4e\t\t%.4e\t\t%.4e\n",nan,p02r,p03r);
    fprintf("Mach, Relative Rotor:\n\t\tM1r\t\t\t\tM2\t\t\t\tM3\n\t\t%.4f\t\t\t%.4f\t\t\t%.4f\n",nan,Mw2,Mw3);
    fprintf("Entropy Change  Values, static  [K]:\n\t\t 1-> 2 \t\t\t\t2 -> 3\n\t\t%.3f\t\t\t\t%.3f\n",nan,Cp*log(t3/t2) - R*log(p3/p2));
    fprintf("Entropy Change  Values, total  Rotor, Relative[K]:\n\t\t 1-> 2 \t\t\t\t2 -> 3\n\t\t%.3f\t\t\t\t%.3f\n",nan,Cp*log(t03r/t02r) - R*log(p03r/p02r));
    
    % Blade Properties - currently table 4, only works when radial
    % resolution is 3
    fprintf("Blade Characteristics\n%s\t%s\n%.4f\t%.4f\n%.4f\t%.4f\n%.4f\t%.4f\n%.4f\t%.4f\n%.4f\t%.4f\n%.4f\t%.4f\n",...
        'Mid Stator','Mid Rotor',cAxialStator_fun(2),cAxialRotor_fun(2),cStator_fun(2),cRotor_fun(2),rad2deg(staggerStator_fun(1)),rad2deg(staggerRotor_fun(1)),...
        rad2deg(staggerStator_fun(2)),rad2deg(staggerRotor_fun(2)),rad2deg(staggerStator_fun(3)),rad2deg(staggerRotor_fun(3)));
    % Meridional View, radii for all positions
    fprintf("Radii for each point:\nCond.1\tS LE\tS TE\tCond.2\tR LE\tR TE\tCond3\n");
    disp(rot90(newRadii));
    fprintf("Area of Annulus for each points:\nCond.1\tS LE\tS TE\tCond.2\tR LE\tR TE\tCond3\n");
    disp(rot90(A_fun));
    fprintf("Corresponding 'x' coordinate of each mid point\n");
    disp(radiusPlotx);
    
    % Geometry / Vel Triangles
    StatorInletVelTriangle =[rad2deg(alpha1_fun); rad2deg(beta1_fun); u1_fun(2,:);v1_fun(1,:);v1_fun(2,:);vecnorm(v1_fun);w1_fun(1,:);w1_fun(2,:);vecnorm(w1_fun)]
    StatorExitRotorInletVelTriangle =abs([rad2deg(alpha2_fun); rad2deg(beta2_fun); u2_fun(2,:);v2_fun(1,:);v2_fun(2,:);vecnorm(v2_fun);w2_fun(1,:);w2_fun(2,:);vecnorm(w2_fun)])
    RotorExitVelTriangle =abs([rad2deg(alpha3_fun); rad2deg(beta3_fun); u3_fun(2,:);v3_fun(1,:);v3_fun(2,:);vecnorm(v3_fun);w3_fun(1,:);w3_fun(2,:);vecnorm(w3_fun)])
%     fprintf("Deflection Angles,alpha [deg] :\n\t\talpha_1\t\t\talpha_2\t\t\talpha_3\n\t\t%.3f\t\t\t%.3f\t\t\t%.3f\n",rad2deg(alpha1),rad2deg(alpha2),rad2deg(alpha3));
%     fprintf("Relative Angles,beta [deg]:\n\t\tbeta_1\t\t\tbeta_2\t\t\tbeta_3\n\t\t-N/A-\t\t\t%.3f\t\t\t%.3f\n",rad2deg(beta2),rad2deg(beta3));
    fprintf("Degree of Reaction at Hub:\n\t\t%.3f%%\n",degOfRxn_fun(1)*100);
    fprintf("Flow Coeff:\n\t\t%.3f\n",phi);
    fprintf("Work Coeff:\n\t\t%.3f\n",lambda);
    fprintf("Mis Span Zweifel Coeff:\n\t\t%.3f\n",zweifel);
    fprintf("Stator Blade Count:\n\t\t%.3f\n",bladeCountStator);
    fprintf("Rotor Blade Count:\n\t\t%.3f\n",bladeCountRotor);
    
fprintf("First Table with design choices:\n");
fprintf("\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n",rpm, zweifel, phi,v2(1)/v1(1)*100,v3(1)/v2(1)*100, rm3/rm2, rm2/rm1, rad2deg(alpha3), rotorAR, rotorTaper, statorAR, statorTaper);

fprintf("Second Table with Major Design Characteristics\n");
fprintf("\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f",powerExtracted/1e6, degOfRxn_fun(1)*100, lambda, maxMach, rad2deg(rotorTurnAng), bladeCountStator, bladeCountRotor);

end

