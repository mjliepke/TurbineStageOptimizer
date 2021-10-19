%% Matthew Liepke, AE 403 Spr 2021
% Script to run through different stage designs and print valid ones.
% Possibly display 3-D visualization of power outputted per stage parameter
clc;clear;close all;
%% Consts & Ranges - [kv1,alpha2, lambda, phi]
% For Desktop Hardware:
%kvRes = 100; % points to test per variable
angleRes = 300;
lambdaRes = 130;
phiRes = 150;
kvRes = 50;
rpmRes = 20;

% For Laptop Hardware / Quick Compiles:
% kvRes = 150; % points to test per variable
% angleRes = 150;
%  otherRes = 50;
mins = [.01,deg2rad(-72), 1.1];
maxes = [1.2,deg2rad(-60),2.1];

% For const when not varying
kv2 = .5051;
rpm = 12750;
phi = .8;
alpha2 = -1.18;
lambda = 2.02;

phi_fun = linspace(.8, 1.1, phiRes);
kv2_fun = linspace(0,1.2,kvRes);
alpha2_fun = linspace(mins(2),maxes(2),angleRes);
lambda_fun  = linspace(mins(3),maxes(3),lambdaRes);
rpm_fun = linspace(11000,13000,rpmRes);

valid_fun = zeros([phiRes,angleRes,lambdaRes]);
power_fun = zeros([phiRes,angleRes,lambdaRes]); %array of all the power variables to plot later (maybe)
maxMach_fun = zeros([phiRes,angleRes,lambdaRes]);
degOfRxnHub_fun = zeros([phiRes,angleRes,lambdaRes]);
statorBlades_fun = zeros([phiRes,angleRes,lambdaRes]);
rotorBlades_fun = zeros([phiRes,angleRes,lambdaRes]);
alpha3_fun = zeros([phiRes,angleRes,lambdaRes]);
rotorTurn_fun = zeros([phiRes,angleRes,lambdaRes]);

tStart = tic;
fprintf("STARTING ITERATIONS\n");
parfor i = 1:phiRes     %Loop through each design parameter combination
    for j=1:angleRes
        for k=1:lambdaRes
                [validStage, powerExtracted, maxMach, degOfRxnHub, statorBladeCount, rotorBladeCount, alpha3, rotorTurnAng] = analyzeStageFunction(kv2,alpha2_fun(j),lambda_fun(k),phi_fun(i), rpm);
               
                if(validStage)
%                    fprintf("VALID STAGE W/ kv1 = %.3f, alpha2 = %.3f, lambda = %.3f, phi = %.3f\n\tPower = %3s, maxMach = %.3f, hubRxn = %.2f\n",...
%                        kv1_fun(i),alpha2_fun(j), lambda_fun(k), phi ,powerExtracted, maxMach, degOfRxnHub*100);
                   valid_fun(i,j,k) = 1;
                   power_fun(i,j,k) = powerExtracted;
                   maxMach_fun(i,j,k) = maxMach;
                   degOfRxnHub_fun(i,j,k) = degOfRxnHub;
                   statorBlades_fun(i,j,k) = statorBladeCount;
                   rotorBlades_fun(i,j,k) = rotorBladeCount;
                   alpha3_fun(i,j,k) = alpha3;
                   rotorTurn_fun(i,j,k) = rotorTurnAng;
                end
                
        end
    end
    fprintf("COMPLETED CASE OF PHI= %f\n",phi_fun(i));
end

tElapsed = toc(tStart)
fprintf("This operation testing %3d cases took %f seconds\n",angleRes*phiRes*lambdaRes, tElapsed);

%% Plots
[is,js,ks] = ind2sub(size(valid_fun),find(valid_fun == 1));
machs = maxMach_fun(sub2ind(size(maxMach_fun),is,js,ks));
powers = power_fun(sub2ind(size(power_fun),is,js,ks));
bladeRatio = statorBlades_fun(sub2ind(size(statorBlades_fun),is,js,ks))./rotorBlades_fun(sub2ind(size(rotorBlades_fun),is,js,ks));
alpha3 = alpha3_fun(sub2ind(size(alpha3_fun),is,js,ks));
rotorTurn = rotorTurn_fun(sub2ind(size(rotorTurn_fun),is,js,ks));

figure('Name','Mach Color Plot');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks), 15, machs,'filled');
xlabel('\Phi');
ylabel('\alpha_2 [rad]');
zlabel('\lambda');
m = colorbar;
ylabel(m,'Max Mach');

figure('Name','Blade Count Ratio');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks),15,  bladeRatio, 'filled');
xlabel('\Phi');
ylabel('\alpha_2 [rad]');
zlabel('\lambda');
m = colorbar;
ylabel(m,'Stator Blade Count / Rotor Blade Count for a Zweifel Coeff of 0.8');

figure('Name','Alpha3');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks),15,  alpha3, 'filled');
xlabel('\Phi');
ylabel('\alpha_2 [rad]');
zlabel('\lambda');
m = colorbar;
ylabel(m,'\alpha_3 [rad]');

figure('Name','Mid Radius Turn Angle');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks),15,  rad2deg(rotorTurn), 'filled');
xlabel('\Phi');
ylabel('\alpha_2 [rad]');
zlabel('\lambda');
m = colorbar;
ylabel(m,'\beta_3+\beta_2 [deg]');

%% Test 2 - Varying the values of phi, lambda and alpha2
valid_fun = zeros([rpm,angleRes,lambdaRes]);

tStart = tic;
fprintf("STARTING ITERATIONS\n");
parfor i = 1:rpmRes     %Loop through each design parameter combination
    for j=1:angleRes
        for k=1:lambdaRes
                [validStage, powerExtracted, maxMach, degOfRxnHub, statorBladeCount, rotorBladeCount, alpha3, rotorTurnAng] =...
                    analyzeStageFunction(kv2,alpha2_fun(j),lambda_fun(k),phi, rpm_fun(i));
               
                if(validStage)
%                    fprintf("VALID STAGE W/ kv1 = %.3f, alpha2 = %.3f, lambda = %.3f, phi = %.3f\n\tPower = %3s, maxMach = %.3f, hubRxn = %.2f\n",...
%                        kv1_fun(i),alpha2_fun(j), lambda_fun(k), phi ,powerExtracted, maxMach, degOfRxnHub*100);
                   valid_fun(i,j,k) = 1;
                   power_fun(i,j,k) = powerExtracted;
                   maxMach_fun(i,j,k) = maxMach;
                   degOfRxnHub_fun(i,j,k) = degOfRxnHub;
                   statorBlades_fun(i,j,k) = statorBladeCount;
                   rotorBlades_fun(i,j,k) = rotorBladeCount;
                   alpha3_fun(i,j,k) = alpha3;
                   rotorTurn_fun(i,j,k) = rotorTurnAng;
                end
                
        end
    end
    fprintf("COMPLETED CASE OF RPM= %f\n",rpm_fun(i));
end

tElapsed = toc(tStart)
fprintf("This operation testing %3d cases took %f seconds\n",angleRes*phiRes*lambdaRes, tElapsed);

%% Plots
[is,js,ks] = ind2sub(size(valid_fun),find(valid_fun == 1));
machs = maxMach_fun(sub2ind(size(maxMach_fun),is,js,ks));
powers = power_fun(sub2ind(size(power_fun),is,js,ks));
bladeRatio = statorBlades_fun(sub2ind(size(statorBlades_fun),is,js,ks))./rotorBlades_fun(sub2ind(size(rotorBlades_fun),is,js,ks));
alpha3 = alpha3_fun(sub2ind(size(alpha3_fun),is,js,ks));
rotorTurn = rotorTurn_fun(sub2ind(size(rotorTurn_fun),is,js,ks));

figure('Name','Mach Color Plot');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks), 15, machs,'filled');
xlabel('rpm');
ylabel('\alpha_2 [rad]');
zlabel('\lambda');
m = colorbar;
ylabel(m,'Max Mach');

figure('Name','Blade Count Ratio');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks),15,  bladeRatio, 'filled');
xlabel('rpm');
ylabel('\alpha_2 [rad]');
zlabel('\lambda');
m = colorbar;
ylabel(m,'Stator Blade Count / Rotor Blade Count for a Zweifel Coeff of 0.8');

figure('Name','Alpha3');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks),15,  alpha3, 'filled');
xlabel('rpm');
ylabel('\alpha_2 [deg]');
zlabel('\lambda');
m = colorbar;
ylabel(m,'\alpha_3 [deg]');

figure('Name','Mid Radius Turn Angle');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks),15,  rotorTurn, 'filled');
xlabel('rpm');
ylabel('\alpha_2 [deg]');
zlabel('\lambda');
m = colorbar;
ylabel(m,'\beta_3+\beta_2 [deg]');


%% Vary RPM, angle, and kv1
valid_fun = zeros([rpmRes,angleRes,kvRes]);

tStart = tic;
fprintf("STARTING ITERATIONS\n");
parfor i = 1:rpmRes     %Loop through each design parameter combination
    for j=1:angleRes
        for k=1:kvRes
                [validStage, powerExtracted, maxMach, degOfRxnHub, statorBladeCount, rotorBladeCount, alpha3, rotorTurnAng] =...
                    analyzeStageFunction(kv2_fun(k),alpha2_fun(j),lambda,phi, rpm_fun(i));
               
                if(validStage)
%                    fprintf("VALID STAGE W/ kv1 = %.3f, alpha2 = %.3f, lambda = %.3f, phi = %.3f\n\tPower = %3s, maxMach = %.3f, hubRxn = %.2f\n",...
%                        kv1_fun(i),alpha2_fun(j), lambda_fun(k), phi ,powerExtracted, maxMach, degOfRxnHub*100);
                   valid_fun(i,j,k) = 1;
                   power_fun(i,j,k) = powerExtracted;
                   maxMach_fun(i,j,k) = maxMach;
                   degOfRxnHub_fun(i,j,k) = degOfRxnHub;
                   statorBlades_fun(i,j,k) = statorBladeCount;
                   rotorBlades_fun(i,j,k) = rotorBladeCount;
                   alpha3_fun(i,j,k) = alpha3;
                   rotorTurn_fun(i,j,k) = rotorTurnAng;
                end
                
        end
    end
    fprintf("COMPLETED CASE OF RPM= %f\n",rpm_fun(i));
end

tElapsed = toc(tStart)
fprintf("This operation testing %3d cases took %f seconds\n",angleRes*phiRes*lambdaRes, tElapsed);

%% Plots
[is,js,ks] = ind2sub(size(valid_fun),find(valid_fun == 1));
machs = maxMach_fun(sub2ind(size(maxMach_fun),is,js,ks));
powers = power_fun(sub2ind(size(power_fun),is,js,ks));
bladeRatio = statorBlades_fun(sub2ind(size(statorBlades_fun),is,js,ks))./rotorBlades_fun(sub2ind(size(rotorBlades_fun),is,js,ks));
alpha3 = alpha3_fun(sub2ind(size(alpha3_fun),is,js,ks));
rotorTurn = rotorTurn_fun(sub2ind(size(rotorTurn_fun),is,js,ks));

figure('Name','Mach Color Plot');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks), 15, machs,'filled');
xlabel('rpm');
ylabel('\alpha_2 [rad]');
zlabel('kv2');
m = colorbar;
ylabel(m,'Max Mach');

figure('Name','Blade Count Ratio');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks),15,  bladeRatio, 'filled');
xlabel('rpm');
ylabel('\alpha_2 [rad]');
zlabel('kv2');
m = colorbar;
ylabel(m,'Stator Blade Count / Rotor Blade Count for a Zweifel Coeff of 0.8');

figure('Name','Alpha3');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks),15,  alpha3, 'filled');
xlabel('rpm');
ylabel('\alpha_2 [deg]');
zlabel('kv2');
m = colorbar;
ylabel(m,'\alpha_3 [deg]');

figure('Name','Mid Radius Turn Angle');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks),15,  rotorTurn, 'filled');
xlabel('rpm');
ylabel('\alpha_2 [deg]');
zlabel('kv2');
m = colorbar;
ylabel(m,'\beta_3+\beta_2 [deg]');