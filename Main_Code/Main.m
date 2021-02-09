% Supplementary Code to the Paper 'A framework for incorporating 3D
% hyperelastic vascular wall models in 1D blood flow simulations' -
% submitted to 'Biomechanics and Modeling in Mechanobiology'.
% Paper Authors - Alberto Coccarelli, Jason M. Carson, Ankush Aggarwal, and Sanjay Pant

% Blood flow Code written by Jason M. Carson (BSc, MSc, PhD), Post-Dotoral Research Fellow at 
% Biomedical Rngineering Research Group, ZCCE, Swansea University
clear
clc

max_sims=2; % number of patient cases run

run('load_carotid_parameters') % Load parameters of carotid artery

%% Selected range for constitutive parameters
H1_range=[0.95;1.25]/10;
Ri_range=[3.5;5.5]/10;
lambzres_range=[0.9;1.1];
lambzE=[1.0];
omega_range=[60;140];

% parameter range
c_range=[min(parameterrange.c);max(parameterrange.c)];
k1_range=[min(parameterrange.k1);max(parameterrange.k1)];
k2_range=[min(parameterrange.k2);max(parameterrange.k2)];
alpha_range=[min(parameterrange.alpha);max(parameterrange.alpha)];
beta_range=[min(parameterrange.beta);max(parameterrange.beta)];

%% Organising and creating a database of constitutive parameters
pv=sobolset(9);
vc=net(pv,max_sims);

H1=H1_range(1)+vc(:,1)*(H1_range(2)-H1_range(1));
Ri=Ri_range(1)+vc(:,2)*(Ri_range(2)-Ri_range(1));
lambzres=lambzres_range(1)+vc(:,3)*(lambzres_range(2)-lambzres_range(1));
omega=omega_range(1)+vc(:,4)*(omega_range(2)-omega_range(1));
c=(c_range(1)+vc(:,5)*(c_range(2)-c_range(1)))*10e4;
k1=(k1_range(1)+vc(:,6)*(k1_range(2)-k1_range(1)))*10e4;
k2=k2_range(1)+vc(:,7)*(k2_range(2)-k2_range(1));
alpha=alpha_range(1)+vc(:,8)*(alpha_range(2)-alpha_range(1));
beta=beta_range(1)+vc(:,9)*(beta_range(2)-beta_range(1));

Set1_IN=[H1,Ri,omega,lambzres,c,k1,k2,alpha,beta]; % organise cases into an array for ease of use 

%% Clear memory
clearvars H1 Ri lambzres omega c k1 k2 alpha beta vc pv c_range k1_range k2_range alpha_range beta_range H1_range Ri_range lambzres_range omega_range parameterrange

%% Loop through cases
for sim_run=1:max_sims
    tic
    % output variables have two dimensions. Rows correspond with the three
    % different constitutive laws. Columns are organised as follows:
% mean pwv using the traditional calculation, mean wavespeed, 
% mean pressure, max pressure, min pressure, mean area, max area, min area,
% difference in pressure between proximal and distal ends, mean compliance,
% maximum compliance, minimum compliance, and mean pwv using the
% intersection method

% in addition the maximum CFL number is recorded for each case (not the
% subdomain collocation scheme implemented is unconditionally stable for
% any CFL number)

    [ Output_Vars,CFL ] = run_three_constitutive_laws([Set1_IN(sim_run,:)],lambzE);
     
    toc
    
    %% Patient number
    disp('Patient Number Finished:=')
    disp(sim_run)
    save(['Results/Results_',num2str(sim_run)],'Output_Vars','CFL' )

end


