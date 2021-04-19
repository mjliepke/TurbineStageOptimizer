%% Matthew Liepke, AE 403 Spr 2021
% Script to run through different stage designs and print valid ones.
% Possibly display 3-D visualization of power outputted per stage parameter
clc;clear;close all;
%% Consts & Ranges - [kv1,alpha2, lambda, phi]
% For Desktop Hardware:
%kvRes = 100; % points to test per variable
angleRes = 300;
lambdaRes = 150;
phiRes = 200;
% For Laptop Hardware / Quick Compiles:
% kvRes = 150; % points to test per variable
% angleRes = 150;
%  otherRes = 50;
mins = [.01,deg2rad(-72), 1.1];
maxes = [1.7,deg2rad(-63),2.1];

phi_fun = linspace(.8, 1.1, phiRes);
kv2 = .5;

%kv1_fun = linspace(mins(1),maxes(1),kvRes);
alpha2_fun = linspace(mins(2),maxes(2),angleRes);
lambda_fun  = linspace(mins(3),maxes(3),lambdaRes);

valid_fun = zeros([phiRes,angleRes,lambdaRes]);
power_fun = zeros([phiRes,angleRes,lambdaRes]); %array of all the power variables to plot later (maybe)
maxMach_fun = zeros([phiRes,angleRes,lambdaRes]);
degOfRxnHub_fun = zeros([phiRes,angleRes,lambdaRes]);
statorBlades_fun = zeros([phiRes,angleRes,lambdaRes]);
rotorBlades_fun = zeros([phiRes,angleRes,lambdaRes]);

tStart = tic;
fprintf("STARTING ITERATIONS\n");
parfor i = 1:phiRes     %Loop through each design parameter combination
    for j=1:angleRes
        for k=1:lambdaRes
                [validStage, powerExtracted, maxMach, degOfRxnHub, statorBladeCount, rotorBladeCount] = analyzeStageFunction(kv2,alpha2_fun(j),lambda_fun(k),phi_fun(i));
               
                if(validStage)
%                    fprintf("VALID STAGE W/ kv1 = %.3f, alpha2 = %.3f, lambda = %.3f, phi = %.3f\n\tPower = %3s, maxMach = %.3f, hubRxn = %.2f\n",...
%                        kv1_fun(i),alpha2_fun(j), lambda_fun(k), phi ,powerExtracted, maxMach, degOfRxnHub*100);
                   valid_fun(i,j,k) = 1;
                   power_fun(i,j,k) = powerExtracted;
                   maxMach_fun(i,j,k) = maxMach;
                   degOfRxnHub_fun(i,j,k) = degOfRxnHub;
                   statorBlades_fun(i,j,k) = statorBladeCount;
                   rotorBlades_fun(i,j,k) = rotorBladeCount;
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


figure('Name','Mach Color Plot');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks), 15, machs,'filled');
xlabel('\Phi');
ylabel('\alpha_2 [rad]');
zlabel('\lambda');
m = colorbar;
ylabel(m,'Max Mach');

figure('Name','Power Production Color Plot');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks),15,  powers, 'filled');
xlabel('\Phi');
ylabel('\alpha_2 [rad]');
zlabel('\lambda');
m = colorbar;
ylabel(m,'Power Produced [W]');

figure('Name','Blade Count Ratio');
scatter3(phi_fun(is),-alpha2_fun(js),lambda_fun(ks),15,  bladeRatio, 'filled');
xlabel('\Phi');
ylabel('\alpha_2 [rad]');
zlabel('\lambda');
m = colorbar;
ylabel(m,'Stator Blade Count / Rotor Blade Count for a Zweifel Coeff of 0.8');

