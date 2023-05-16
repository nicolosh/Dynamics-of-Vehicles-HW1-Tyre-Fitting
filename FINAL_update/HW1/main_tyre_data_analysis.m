%% Initialisation
clc
clearvars 
close all   

% Set LaTeX as default interpreter for axis labels, ticks and legends
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultAxesFontSize',  16)
set(0,'DefaultLegendFontSize',16)

addpath('tyre_lib/')


to_rad = pi/180;
to_deg = 180/pi;

%% LATERAL
%% Select tyre dataset
%dataset path
data_set_path = 'dataset/';

% dataset selection and loading
data_set = 'Hoosier_B1464run23'; % pure lateral forces
%data_set = 'Hoosier_B1464run30';  % braking/traction (pure log. force) + combined

% tyre geometric data:
% Hoosier	18.0x6.0-10
% 18 diameter in inches
% 6.0 section width in inches
% tread width in inches
diameter = 18*2.56; %
Fz0 = 1120;   % [N] nominal load is given
R0  = diameter/2/100; % [m] get from nominal load R0 (m) *** TO BE CHANGED ***

fprintf('Loading dataset ...')
switch data_set
    case 'Hoosier_B1464run23'
  load ([data_set_path, 'Hoosier_B1464run23.mat']); % pure lateral
  cut_start = 31300;
  cut_end   = 54500;
    case 'Hoosier_B1464run30'
  load ([data_set_path, 'Hoosier_B1464run30.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 37643;
  otherwise 
  error('Not found dataset: `%s`\n', data_set) ;
end

% select dataset portion
smpl_range = cut_start:cut_end;

fprintf('completed!\n')
%% Plot raw data

figure
tiledlayout(6,1)

ax_list(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list,'x')

%% Select some specific data
% Cut crappy data and select only 12 psi data
clear idx;
vec_samples = 1:1:length(smpl_range);
tyre_data = table(); % create empty table
% store raw data in table
tyre_data.SL =  SL(smpl_range);
tyre_data.SA =  SA(smpl_range)*to_rad;
tyre_data.FZ = -FZ(smpl_range);  % 0.453592  lb/kg
tyre_data.FX =  FX(smpl_range);
tyre_data.FY =  FY(smpl_range);
tyre_data.MZ =  MZ(smpl_range);
tyre_data.IA =  IA(smpl_range)*to_rad;
% Cutting out parts where SA is accumulating (around 0 and -0.068 rads)
% SA_tol = 0.07*to_rad;
% idx.cuttingzeros = 0.0*to_rad-SA_tol < tyre_data.SA & tyre_data.SA < 0.0*to_rad+SA_tol;
% idx.cutting0068 =  -0.068-SA_tol < tyre_data.SA & tyre_data.SA < -0.068+SA_tol;
% idx.cuttingcrap = ~idx.cuttingzeros & ~idx.cutting0068;
% tyre_data = tyre_data( idx.cuttingcrap, : );
%tyre_data.SA = -tyre_data.SA;
% Extract points at constant inclination angle (camber)
GAMMA_tol = 0.05*to_rad;
idx.GAMMA_0 = 0.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 0.0*to_rad+GAMMA_tol;
idx.GAMMA_1 = 1.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 1.0*to_rad+GAMMA_tol;
idx.GAMMA_2 = 2.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 2.0*to_rad+GAMMA_tol;
idx.GAMMA_3 = 3.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 3.0*to_rad+GAMMA_tol;
idx.GAMMA_4 = 4.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 4.0*to_rad+GAMMA_tol;
idx.GAMMA_5 = 5.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 5.0*to_rad+GAMMA_tol;
for i = 310:565:length(tyre_data.IA)
    k=i;
    for j = k:1:k+320
        idx.GAMMA_VAR(j) = 1&1;
        i = j;
        if j == length(tyre_data.IA)
            break
        end
    end
end
idx.GAMMA_VAR = idx.GAMMA_VAR.';
GAMMA_0  = tyre_data( idx.GAMMA_0, : );
GAMMA_1  = tyre_data( idx.GAMMA_1, : );
GAMMA_2  = tyre_data( idx.GAMMA_2, : );
GAMMA_3  = tyre_data( idx.GAMMA_3, : );
GAMMA_4  = tyre_data( idx.GAMMA_4, : );
GAMMA_5  = tyre_data( idx.GAMMA_5, : );
GAMMA_VAR = tyre_data( idx.GAMMA_VAR, : );
% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

FZ_tol = 100;
idx.FZ_220  = 220-FZ_tol < tyre_data.FZ & tyre_data.FZ < 220+FZ_tol;
idx.FZ_440  = 440-FZ_tol < tyre_data.FZ & tyre_data.FZ < 440+FZ_tol;
idx.FZ_700  = 700-FZ_tol < tyre_data.FZ & tyre_data.FZ < 700+FZ_tol;
idx.FZ_900  = 900-FZ_tol < tyre_data.FZ & tyre_data.FZ < 900+FZ_tol;
idx.FZ_1120 = 1120-FZ_tol < tyre_data.FZ & tyre_data.FZ < 1120+FZ_tol;
FZ_220  = tyre_data( idx.FZ_220, : );
FZ_440  = tyre_data( idx.FZ_440, : );
FZ_700  = tyre_data( idx.FZ_700, : );
FZ_900  = tyre_data( idx.FZ_900, : );
FZ_1120 = tyre_data( idx.FZ_1120, : );


figure()
tiledlayout(3,1)

ax_list(1) = nexttile;
plot(tyre_data.IA*to_deg)
hold on
plot(vec_samples(idx.GAMMA_0),GAMMA_0.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_1),GAMMA_1.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_2),GAMMA_2.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_3),GAMMA_3.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_4),GAMMA_4.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_5),GAMMA_5.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(2) = nexttile;
plot(tyre_data.FZ)
hold on
plot(vec_samples(idx.FZ_220),FZ_220.FZ,'.');
plot(vec_samples(idx.FZ_440),FZ_440.FZ,'.');
plot(vec_samples(idx.FZ_700),FZ_700.FZ,'.');
plot(vec_samples(idx.FZ_900),FZ_900.FZ,'.');
plot(vec_samples(idx.FZ_1120),FZ_1120.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')


ax_list(3) = nexttile;
plot(tyre_data.SA)
hold on
title('Slide slip')
xlabel('Samples [-]')
ylabel('[rad]')

%% Intersect tables to obtain specific sub-datasets
[TData1, ~] = intersect_table_data( GAMMA_0, FZ_1120 );

%% initialise tyre data
tyre_coeffs = initialise_tyre_data(R0, Fz0);

%% Fitting with Fz=Fz_nom= 1120N and camber=0 VX= 10
FZ0 = mean(TData1.FZ);
ALPHA_vec = TData1.SA;
FY_vec = TData1.FY;

% plot raw data 
figure('Name','FY0 guess')
plot(TData1.SA,TData1.FY,'.')

% Guess values for parameters to be optimised
%    [pCy1 pDy1 pEy1 pHy1  pKy1  pKy2  pVy1] 
% P0 = [1.5, 2.5, 0.3, 0.1, 100, 1, 0.1];
P0 = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
%    [pCy1 pDy1 pEy1 pEy3 pHy1 pKy1 pKy2  pVy1]
lb = [  1,-1000,  -1000,  -1000, -1000, -1000, -1000];
ub = [  2, 1000,   1000,   1000,  1000,  1000, 1000];
% resid_pure_Fy returns the residual, so minimize the residual varying P. It
% is an unconstrained minimization problem 
options = optimset('MaxFunEvals',3.0e+10);
[P_fz_nom,fval,exitflag] = fmincon(@(P)resid_pure_Fy(P,FY_vec, ALPHA_vec,0,FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub,[],options);


% Update tyre_coeffs with new optimal values                             
tyre_coeffs.pCy1 = P_fz_nom(1) ; 
tyre_coeffs.pDy1 = P_fz_nom(2) ;  
tyre_coeffs.pEy1 = P_fz_nom(3) ;
tyre_coeffs.pHy1 = P_fz_nom(4) ;
tyre_coeffs.pKy1 = P_fz_nom(5) ; 
tyre_coeffs.pKy2 = P_fz_nom(6) ;
tyre_coeffs.pVy1 = P_fz_nom(7) ;

R2.FY_FZnom = 1-resid_pure_Fy(P_fz_nom, FY_vec, ALPHA_vec,0,FZ0,tyre_coeffs);
RMSE.FY_FZnom = sqrt((resid_pure_Fy(P_fz_nom, FY_vec, ALPHA_vec,0,FZ0,tyre_coeffs)*sum(FY_vec.^2))/length(FY_vec));


% Plot fitted curve and raw data
SA_vec = -0.3:0.001:0.3;
FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                              FZ0.*ones(size(SA_vec)),tyre_coeffs);
figure('Name','Fy0(Fz0)')
plot(TData1.SA,TData1.FY,'o')
hold on
plot(SA_vec,FY0_fz_nom_vec,'-','LineWidth',2)
xlabel('$\alpha$ [rad]')
ylabel('$F_{y0}$ [N]')
title('Pure lateral force - Raw data vs fitted curve - FZ = 1120 [N] $\gamma = 0$ [deg]')
legend('Raw data', 'Fitted curve')

%% Fit coefficient with variable load
%[TDataFzVar, ~] = intersect_table_data( GAMMA_0, FZ_1120);
ALPHA_vec = GAMMA_0.SA;
FY_vec    = GAMMA_0.FY;
FZ_vec    = GAMMA_0.FZ;

% Plot raw data 
figure()
plot(GAMMA_0.SA,GAMMA_0.FY,'.')

% Guess values for parameters to be optimised
%    [ pDy2   pEy2   pHy2   pVy2 ] 
P0 = [ 0.1,   0.1,   0.1,    0.1 ]; 
lb = [-1000,   0,   -1000, -1000 ];
ub = [ 1000,   1,    1000,  1000 ];
% minimization
[P_dfz,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varFz(P,FY_vec, ALPHA_vec,0,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre_coeffs with new optimal values                             
tyre_coeffs.pDy2 = P_dfz(1) ; 
tyre_coeffs.pEy2 = P_dfz(2) ;  
tyre_coeffs.pHy2 = P_dfz(3) ;
tyre_coeffs.pVy2 = P_dfz(4) ;

% R2 & RMSE
R2.FY_FZvar = 1-resid_pure_Fy_varFz(P_dfz, FY_vec, ALPHA_vec,0,FZ_vec,tyre_coeffs);
RMSE.FY_FZvar = sqrt((resid_pure_Fy_varFz(P_dfz, FY_vec, ALPHA_vec,0,FZ_vec,tyre_coeffs)*sum(FY_vec.^2))/length(FY_vec));

% Compute residuals with new coeffs
SA_vec = -0.3:0.001:0.3;
res_FY0_dfz_vec = resid_pure_Fy_varFz(P_dfz,FY_vec,SA_vec,0 , FZ_vec,tyre_coeffs);

% Plot curves
dataplot1 = intersect_table_data(GAMMA_0, FZ_220);
dataplot2 = intersect_table_data(GAMMA_0, FZ_440);
dataplot3 = intersect_table_data(GAMMA_0, FZ_700);
dataplot4 = intersect_table_data(GAMMA_0, FZ_900);
dataplot5 = intersect_table_data(GAMMA_0, FZ_1120);
tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));
FY0_fz_var_vec1 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec2 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_440.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec3 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec4 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec5 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);
figure('Name','Fy0(Fz)')
plot(GAMMA_0.SA,GAMMA_0.FY,'o')
hold on
plot(SA_vec,FY0_fz_var_vec1,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec2,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec3,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec4,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec5,'-','LineWidth',2)
legend('Raw data','FZ = 220', 'FZ = 440', 'FZ = 700', 'FZ = 900', 'FZ = 1120')
xlabel('$\alpha$ [rad]')
ylabel('$F_{y0}$ [N]')
title('Pure lateral force - Raw data vs fitted curves - FZ = var [N] $\gamma = 0$ [deg]')

%% Fit coefficient with variable camber
% Cutting SA=0 accumulation
SA_tol = 0.09*to_rad;
idx.cuttingzeros = 0.0*to_rad-SA_tol < FZ_1120.SA & FZ_1120.SA < 0.0*to_rad+SA_tol;
FZ_1120_nozero = FZ_1120( ~idx.cuttingzeros, : );

zeros_vec = zeros(size(FZ_1120_nozero));
ones_vec  = ones(size(FZ_1120_nozero));
ALPHA_vec = FZ_1120_nozero.SA;
GAMMA_vec = FZ_1120_nozero.IA; 
FY_vec    = FZ_1120_nozero.FY;

% Plot raw data
figure()
plot(ALPHA_vec,FY_vec);

% Guess values for parameters to be optimised
% [pDy3 pEy3 pEy4 pHy3 pKy3 pVy3 ]
P0 = [0.1, 0.5, 0.1, 0.1, 0.1, 5]; 
lb = [];
ub = [];
 
[P_varGamma1,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varGamma(P,FY_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre_coeffs with new optimal values    
tyre_coeffs.pDy3 = P_varGamma1(1) ; 
tyre_coeffs.pEy3 = P_varGamma1(2) ;
tyre_coeffs.pEy4 = P_varGamma1(3) ;
tyre_coeffs.pHy3 = P_varGamma1(4) ;
tyre_coeffs.pKy3 = P_varGamma1(5) ;
tyre_coeffs.pVy3 = P_varGamma1(6) ;

% R2 & RMSE
R2.FY_GAMMAvar = 1-resid_pure_Fy_varGamma(P_varGamma1, FY_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs.FZ0,tyre_coeffs);
RMSE.FY_GAMMAvar = sqrt((resid_pure_Fy_varGamma(P_varGamma1, FY_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs.FZ0,tyre_coeffs)*sum(FY_vec.^2))/length(FY_vec));

% Plots
FY0_varGamma_vec0 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_0.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
FY0_varGamma_vec1 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_1.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
FY0_varGamma_vec2 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_2.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
FY0_varGamma_vec3 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_3.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
FY0_varGamma_vec4 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_4.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

figure()
hold on
plot(intersect_table_data(GAMMA_0, FZ_220).SA,intersect_table_data(GAMMA_0, FZ_220).FY,'o','Color','r')
plot(intersect_table_data(GAMMA_4, FZ_220).SA,intersect_table_data(GAMMA_4, FZ_220).FY,'o','Color','b')
legend('$\gamma$ = 0','$\gamma$ = 4')

figure('Name','Fy0 vs Gamma')
hold on
plot(ALPHA_vec,FZ_1120_nozero.FY,'o')
plot(SA_vec,FY0_varGamma_vec0,'-')
plot(SA_vec,FY0_varGamma_vec1,'-')
plot(SA_vec,FY0_varGamma_vec2,'-')
plot(SA_vec,FY0_varGamma_vec3,'-')
plot(SA_vec,FY0_varGamma_vec4,'-')
xlabel('$\alpha$ [rad]')
ylabel('$F_{y0}$ [N]')
legend('Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$')
title('Pure lateral force - Raw data vs fitted curves - FZ = 1120 [N] $\gamma = var$ [deg]')
figure()
hold on
tiledlayout(3,2)

nexttile
hold on
plot(intersect_table_data(GAMMA_0, FZ_1120_nozero).SA,intersect_table_data(GAMMA_0, FZ_1120_nozero).FY,'o')
plot(SA_vec,FY0_varGamma_vec0,'-')
title('$\gamma = 0$ Raw vs fitted')

nexttile
hold on
plot(intersect_table_data(GAMMA_1, FZ_1120_nozero).SA,intersect_table_data(GAMMA_1, FZ_1120_nozero).FY,'o')
plot(SA_vec,FY0_varGamma_vec1,'-')
title('$\gamma = 1$ Raw vs fitted')

nexttile
hold on
plot(intersect_table_data(GAMMA_2, FZ_1120_nozero).SA,intersect_table_data(GAMMA_2, FZ_1120_nozero).FY,'o')
plot(SA_vec,FY0_varGamma_vec2,'-')
title('$\gamma = 2$ Raw vs fitted')

nexttile
hold on
plot(intersect_table_data(GAMMA_3, FZ_1120_nozero).SA,intersect_table_data(GAMMA_3, FZ_1120_nozero).FY,'o')
plot(SA_vec,FY0_varGamma_vec3,'-')
title('$\gamma = 3$ Raw vs fitted')

nexttile
hold on
plot(intersect_table_data(GAMMA_4, FZ_1120_nozero).SA,intersect_table_data(GAMMA_4, FZ_1120_nozero).FY,'o')
plot(SA_vec,FY0_varGamma_vec4,'-')
title('$\gamma = 4$ Raw vs fitted')


% Cornering stiffness
Calfa_vec1 = MF96_CorneringStiffness2(tmp_zeros,SA_vec,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec2 = MF96_CorneringStiffness2(tmp_zeros,SA_vec,tmp_zeros, mean(FZ_440.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffness2(tmp_zeros,SA_vec,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec4 = MF96_CorneringStiffness2(tmp_zeros,SA_vec,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec5 = MF96_CorneringStiffness2(tmp_zeros,SA_vec,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

figure('Name','C_alpha')
hold on
plot(SA_vec,Calfa_vec1,'-','LineWidth',2)
plot(SA_vec,Calfa_vec2,'-','LineWidth',2)
plot(SA_vec,Calfa_vec3,'-','LineWidth',2)
plot(SA_vec,Calfa_vec4,'-','LineWidth',2)
plot(SA_vec,Calfa_vec5,'-','LineWidth',2)
xlabel('$\alpha$ [rad]')
ylabel('$C_{F_\alpha}$ [N/rad]')
legend({'$Fz_{220}$','$Fz_{440}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})
title('Cornering stiffness')

%% FZ = var, gamma = var

ALPHA_vec = tyre_data.SA;
FY_vec = tyre_data.FY;
FZ_vec = tyre_data.FZ;
GAMMA_vec = tyre_data.IA;

% Guess values for parameters to be optimised
%    [ pVy4 ]
P0 = 0;
lb = [];
ub = [];

[P_FY0_dfz_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varFz_varGamma(P,FY_vec, ALPHA_vec,GAMMA_vec,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre_coeffs with new optimal values 
tyre_coeffs.pVy4 = P_FY0_dfz_varGamma(1);

% R2 & RMSE
R2.FY_FZvar_GAMMAvar = 1-resid_pure_Fy_varFz_varGamma(P_FY0_dfz_varGamma, FY_vec, ALPHA_vec,GAMMA_vec,FZ_vec,tyre_coeffs);
RMSE.FY_FZvar_GAMMAvar = sqrt((resid_pure_Fy_varFz_varGamma(P_FY0_dfz_varGamma, FY_vec, ALPHA_vec,GAMMA_vec,FZ_vec,tyre_coeffs)*sum(FY_vec.^2))/length(FY_vec));

% Plots
SA_vec = -0.3:0.001:0.3;
FY0_dfz_varGamma_vec1 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_0.IA,mean(FZ_220.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec2 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_1.IA,mean(FZ_220.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec3 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_2.IA,mean(FZ_220.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec4 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_3.IA,mean(FZ_220.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec5 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_4.IA,mean(FZ_220.FZ)*ones(size(SA_vec)),tyre_coeffs);

FY0_dfz_varGamma_vec6 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_0.IA,mean(FZ_440.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec7 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_1.IA,mean(FZ_440.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec8 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_2.IA,mean(FZ_440.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec9 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_3.IA,mean(FZ_440.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec10 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_4.IA,mean(FZ_440.FZ)*ones(size(SA_vec)),tyre_coeffs);

FY0_dfz_varGamma_vec11 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_0.IA,mean(FZ_700.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec12 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_1.IA,mean(FZ_700.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec13 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_2.IA,mean(FZ_700.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec14 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_3.IA,mean(FZ_700.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec15 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_4.IA,mean(FZ_700.FZ)*ones(size(SA_vec)),tyre_coeffs);

FY0_dfz_varGamma_vec16 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_0.IA,mean(FZ_900.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec17 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_1.IA,mean(FZ_900.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec18 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_2.IA,mean(FZ_900.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec19 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_3.IA,mean(FZ_900.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec20 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_4.IA,mean(FZ_900.FZ)*ones(size(SA_vec)),tyre_coeffs);

FY0_dfz_varGamma_vec21 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_0.IA,mean(FZ_1120.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec22 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_1.IA,mean(FZ_1120.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec23 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_2.IA,mean(FZ_1120.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec24 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_3.IA,mean(FZ_1120.FZ)*ones(size(SA_vec)),tyre_coeffs);
FY0_dfz_varGamma_vec25 = MF96_FY0_vec(zeros(size(SA_vec)),SA_vec,GAMMA_4.IA,mean(FZ_1120.FZ)*ones(size(SA_vec)),tyre_coeffs);

figure()
title('Self aligning moment - Raw data vs fitted curve - $\gamma = var$, $FZ = var$')
hold on
plot(ALPHA_vec,FY_vec)
plot(SA_vec,FY0_dfz_varGamma_vec1)
plot(SA_vec,FY0_dfz_varGamma_vec2)
plot(SA_vec,FY0_dfz_varGamma_vec3)
plot(SA_vec,FY0_dfz_varGamma_vec4)
plot(SA_vec,FY0_dfz_varGamma_vec5)
plot(SA_vec,FY0_dfz_varGamma_vec6)
plot(SA_vec,FY0_dfz_varGamma_vec7)
plot(SA_vec,FY0_dfz_varGamma_vec8)
plot(SA_vec,FY0_dfz_varGamma_vec9)
plot(SA_vec,FY0_dfz_varGamma_vec10)
plot(SA_vec,FY0_dfz_varGamma_vec11)
plot(SA_vec,FY0_dfz_varGamma_vec12)
plot(SA_vec,FY0_dfz_varGamma_vec13)
plot(SA_vec,FY0_dfz_varGamma_vec14)
plot(SA_vec,FY0_dfz_varGamma_vec15)
plot(SA_vec,FY0_dfz_varGamma_vec16)
plot(SA_vec,FY0_dfz_varGamma_vec17)
plot(SA_vec,FY0_dfz_varGamma_vec18)
plot(SA_vec,FY0_dfz_varGamma_vec19)
plot(SA_vec,FY0_dfz_varGamma_vec20)
plot(SA_vec,FY0_dfz_varGamma_vec21)
plot(SA_vec,FY0_dfz_varGamma_vec22)
plot(SA_vec,FY0_dfz_varGamma_vec23)
plot(SA_vec,FY0_dfz_varGamma_vec24)
plot(SA_vec,FY0_dfz_varGamma_vec25)


figure()
hold on
tiledlayout(2,3)
title('Pure lateral force - Raw data vs fitted curves - $\gamma = var$, $FZ = var$')
nexttile
hold on
plot(FZ_220.SA,FZ_220.FY,'o')
plot(SA_vec,FY0_dfz_varGamma_vec1)
plot(SA_vec,FY0_dfz_varGamma_vec2)
plot(SA_vec,FY0_dfz_varGamma_vec3)
plot(SA_vec,FY0_dfz_varGamma_vec4)
plot(SA_vec,FY0_dfz_varGamma_vec5)
title('FZ = 220')
legend('Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$')
xlabel('$\alpha$ [rad]')
ylabel('$F_{Y0} [N]$')
nexttile
hold on
plot(FZ_440.SA,FZ_440.FY,'o')
plot(SA_vec,FY0_dfz_varGamma_vec6)
plot(SA_vec,FY0_dfz_varGamma_vec7)
plot(SA_vec,FY0_dfz_varGamma_vec8)
plot(SA_vec,FY0_dfz_varGamma_vec9)
plot(SA_vec,FY0_dfz_varGamma_vec10)
title('FZ = 440')
legend('Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$')
xlabel('$\alpha$ [rad]')
ylabel('$F_{Y0} [N]$')
nexttile
hold on
plot(FZ_700.SA,FZ_700.FY,'o')
plot(SA_vec,FY0_dfz_varGamma_vec11)
plot(SA_vec,FY0_dfz_varGamma_vec12)
plot(SA_vec,FY0_dfz_varGamma_vec13)
plot(SA_vec,FY0_dfz_varGamma_vec14)
plot(SA_vec,FY0_dfz_varGamma_vec15)
title('FZ = 700')
legend('Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$')
xlabel('$\alpha$ [rad]')
ylabel('$F_{Y0} [N]$')
nexttile
hold on
plot(FZ_900.SA,FZ_900.FY,'o')
plot(SA_vec,FY0_dfz_varGamma_vec16)
plot(SA_vec,FY0_dfz_varGamma_vec17)
plot(SA_vec,FY0_dfz_varGamma_vec18)
plot(SA_vec,FY0_dfz_varGamma_vec19)
plot(SA_vec,FY0_dfz_varGamma_vec20)
title('FZ = 900')
legend('Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$')
xlabel('$\alpha$ [rad]')
ylabel('$F_{Y0} [N]$')
nexttile
hold on
plot(FZ_1120.SA,FZ_1120.FY,'o')
plot(SA_vec,FY0_dfz_varGamma_vec21)
plot(SA_vec,FY0_dfz_varGamma_vec22)
plot(SA_vec,FY0_dfz_varGamma_vec23)
plot(SA_vec,FY0_dfz_varGamma_vec24)
plot(SA_vec,FY0_dfz_varGamma_vec25)
title('FZ = 1120')
legend('Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$')
xlabel('$\alpha$ [rad]')
ylabel('$F_{Y0} [N]$')

%% Self Aligning Moment
%% FZ = 1120, gamma = 0

FZ0 = mean(TData1.FZ);
ALPHA_vec = TData1.SA;
MZ_vec = TData1.MZ;

% Guess values for parameters to be optimised
%    [ qHz1 qBz1 qCz1 qDz1 qEz1 qEz4 qBz9 qBz10 qDz6 ] 

P0 = [  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,  0.1, 0.1];
lb = [];
ub = [];

[P_MZ0,fval,exitflag] = fmincon(@(P)resid_Mz(P,MZ_vec, ALPHA_vec,0,FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre_coeffs with new optimal values 
tyre_coeffs.qHz1 = P_MZ0(1);
tyre_coeffs.qBz1 = P_MZ0(2);
tyre_coeffs.qCz1 = P_MZ0(3);
tyre_coeffs.qDz1 = P_MZ0(4);
tyre_coeffs.qEz1 = P_MZ0(5);
tyre_coeffs.qEz4 = P_MZ0(6);
tyre_coeffs.qBz9 = P_MZ0(7);
tyre_coeffs.qBz10 = P_MZ0(8);
tyre_coeffs.qDz6 = P_MZ0(9);

% R2 & RMSE
R2.MZ_FZnom = 1-resid_Mz(P_MZ0, MZ_vec, ALPHA_vec,0,FZ0,tyre_coeffs);
RMSE.MZ_FZnom = sqrt((resid_Mz(P_MZ0, MZ_vec, ALPHA_vec,0,FZ0,tyre_coeffs)*sum(MZ_vec.^2))/length(MZ_vec));

% Plots
MZ0_vec = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,zeros(size(SA_vec)),1120*ones(size(SA_vec)),tyre_coeffs);
figure()
title('Self aligning moment - Raw data vs fitted curve - $\gamma = 0$, $FZ = 1120$')
hold on
plot(ALPHA_vec,MZ_vec,'o')
plot(SA_vec,MZ0_vec,'LineWidth',2)
xlabel('$\alpha$ [rad]')
ylabel('$M_{z} [N \cdot m]$')
%% FZ = var, gamma = 0

ALPHA_vec = GAMMA_0.SA;
MZ_vec = GAMMA_0.MZ;
FZ_vec = GAMMA_0.FZ;
% Guess values for parameters to be optimised
%    [ qHz2 qBz2 qBz3 qDz2 qEz2 qEz3 qDz7 ] 
P0 = [  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ];
lb = [];
ub = [];

[P_MZ0_dfz,fval,exitflag] = fmincon(@(P)resid_Mz_dfz(P,MZ_vec, ALPHA_vec,0,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre_coeffs with new optimal values 
tyre_coeffs.qHz2 = P_MZ0_dfz(1);
tyre_coeffs.qBz2 = P_MZ0_dfz(2);
tyre_coeffs.qBz3 = P_MZ0_dfz(3);
tyre_coeffs.qDz2 = P_MZ0_dfz(4);
tyre_coeffs.qEz2 = P_MZ0_dfz(5);
tyre_coeffs.qEz3 = P_MZ0_dfz(6);
tyre_coeffs.qDz7 = P_MZ0_dfz(7);

% R2 & RMSE
R2.MZ_FZvar = 1-resid_Mz_dfz(P_MZ0_dfz, MZ_vec, ALPHA_vec,0,FZ_vec,tyre_coeffs);
RMSE.MZ_FZvar = sqrt((resid_Mz_dfz(P_MZ0_dfz, MZ_vec, ALPHA_vec,0,FZ_vec,tyre_coeffs)*sum(MZ_vec.^2))/length(MZ_vec));

% Plots
MZ0_dfz_vec1 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,zeros(size(SA_vec)),mean(FZ_220.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_vec2 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,zeros(size(SA_vec)),mean(FZ_440.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_vec3 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,zeros(size(SA_vec)),mean(FZ_700.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_vec4 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,zeros(size(SA_vec)),mean(FZ_900.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_vec5 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,zeros(size(SA_vec)),mean(FZ_1120.FZ)*ones(size(SA_vec)),tyre_coeffs);

figure()
title('Self aligning moment - Raw data vs fitted curves - $\gamma = 0$, $FZ = var$')
hold on
plot(ALPHA_vec,GAMMA_0.MZ,'o')
plot(SA_vec,MZ0_dfz_vec1,'LineWidth',2)
plot(SA_vec,MZ0_dfz_vec2,'LineWidth',2)
plot(SA_vec,MZ0_dfz_vec3,'LineWidth',2)
plot(SA_vec,MZ0_dfz_vec4,'LineWidth',2)
plot(SA_vec,MZ0_dfz_vec5,'LineWidth',2)
legend('Raw data','FZ = 220','FZ = 440','FZ = 700','FZ = 900','FZ = 1120')
title('Self aligning moment - Raw data vs fitted curves - $\gamma = 0$, $FZ = var$')

%% FZ = 1120, gamma = var

ALPHA_vec = FZ_1120.SA;
MZ_vec = FZ_1120.MZ;
FZ_mean = mean(FZ_1120.FZ);
GAMMA_vec = FZ_1120.IA;

% Guess values for parameters to be optimised
%    [ qHz3 qBz4 qBz5 qDz3 qDz4 qEz5 qDz8 ] 
P0 = [  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ];
lb = [];
ub = [];

[P_MZ0_varGamma,fval,exitflag] = fmincon(@(P)resid_Mz_varGamma(P,MZ_vec, ALPHA_vec,GAMMA_vec,FZ_mean, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre_coeffs with new optimal values 
tyre_coeffs.qHz3 = P_MZ0_varGamma(1);
tyre_coeffs.qBz4 = P_MZ0_varGamma(2);
tyre_coeffs.qBz5 = P_MZ0_varGamma(3);
tyre_coeffs.qDz3 = P_MZ0_varGamma(4);
tyre_coeffs.qDz4 = P_MZ0_varGamma(5);
tyre_coeffs.qEz5 = P_MZ0_varGamma(6);
tyre_coeffs.qDz8 = P_MZ0_varGamma(7);

% R2 & RMSE
R2.MZ_GAMMAvar = 1-resid_Mz_varGamma(P_MZ0_varGamma, MZ_vec, ALPHA_vec,GAMMA_vec,FZ_mean,tyre_coeffs);
RMSE.MZ_GAMMAvar = sqrt((resid_Mz_varGamma(P_MZ0_varGamma, MZ_vec, ALPHA_vec,GAMMA_vec,FZ_mean,tyre_coeffs)*sum(MZ_vec.^2))/length(MZ_vec));

% Plots
MZ0_varGamma_vec1 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_0.IA,FZ_mean*ones(size(SA_vec)),tyre_coeffs);
MZ0_varGamma_vec2 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_1.IA,FZ_mean*ones(size(SA_vec)),tyre_coeffs);
MZ0_varGamma_vec3 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_2.IA,FZ_mean*ones(size(SA_vec)),tyre_coeffs);
MZ0_varGamma_vec4 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_3.IA,FZ_mean*ones(size(SA_vec)),tyre_coeffs);
MZ0_varGamma_vec5 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_4.IA,FZ_mean*ones(size(SA_vec)),tyre_coeffs);

figure()
title('Self aligning moment - Raw data vs fitted curves - $\gamma = var$, $FZ = 1120$')
hold on
plot(ALPHA_vec,MZ_vec,'o')
plot(SA_vec,MZ0_varGamma_vec1,'LineWidth',2)
plot(SA_vec,MZ0_varGamma_vec2,'LineWidth',2)
plot(SA_vec,MZ0_varGamma_vec3,'LineWidth',2)
plot(SA_vec,MZ0_varGamma_vec4,'LineWidth',2)
plot(SA_vec,MZ0_varGamma_vec5,'LineWidth',2)
legend('Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$')

%% FZ = var, gamma = var

ALPHA_vec = tyre_data.SA;
MZ_vec = tyre_data.MZ;
FZ_vec = tyre_data.FZ;
GAMMA_vec = tyre_data.IA;

% Guess values for parameters to be optimised
%    [ qHz4 qDz9 ] 
P0 = [ 0.5,    1 ];
lb = [ 0.4,  0.9 ];
ub = [ 0.6,  1.1 ];

[P_MZ0_dfz_varGamma,fval,exitflag] = fmincon(@(P)resid_Mz_dfz_varGamma(P,MZ_vec, ALPHA_vec,GAMMA_vec,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre_coeffs with new optimal values 
tyre_coeffs.qHz4 = P_MZ0_dfz_varGamma(1);
tyre_coeffs.qDz9 = P_MZ0_dfz_varGamma(2);

% R2 & RMSE
R2.MZ_FZvar_GAMMAvar = 1-resid_Mz_dfz_varGamma(P_MZ0_dfz_varGamma, MZ_vec, ALPHA_vec,GAMMA_vec,FZ_vec,tyre_coeffs);
RMSE.MZ_FZvar_GAMMAvar = sqrt((resid_Mz_dfz_varGamma(P_MZ0_dfz_varGamma, MZ_vec, ALPHA_vec,GAMMA_vec,FZ_vec,tyre_coeffs)*sum(MZ_vec.^2))/length(MZ_vec));

% Plots
MZ0_dfz_varGamma_vec1 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_0.IA,mean(FZ_220.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec2 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_1.IA,mean(FZ_220.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec3 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_2.IA,mean(FZ_220.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec4 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_3.IA,mean(FZ_220.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec5 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_4.IA,mean(FZ_220.FZ)*ones(size(SA_vec)),tyre_coeffs);

MZ0_dfz_varGamma_vec6 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_0.IA,mean(FZ_440.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec7 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_1.IA,mean(FZ_440.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec8 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_2.IA,mean(FZ_440.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec9 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_3.IA,mean(FZ_440.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec10 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_4.IA,mean(FZ_440.FZ)*ones(size(SA_vec)),tyre_coeffs);

MZ0_dfz_varGamma_vec11 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_0.IA,mean(FZ_700.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec12 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_1.IA,mean(FZ_700.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec13 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_2.IA,mean(FZ_700.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec14 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_3.IA,mean(FZ_700.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec15 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_4.IA,mean(FZ_700.FZ)*ones(size(SA_vec)),tyre_coeffs);

MZ0_dfz_varGamma_vec16 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_0.IA,mean(FZ_900.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec17 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_1.IA,mean(FZ_900.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec18 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_2.IA,mean(FZ_900.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec19 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_3.IA,mean(FZ_900.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec20 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_4.IA,mean(FZ_900.FZ)*ones(size(SA_vec)),tyre_coeffs);

MZ0_dfz_varGamma_vec21 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_0.IA,mean(FZ_1120.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec22 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_1.IA,mean(FZ_1120.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec23 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_2.IA,mean(FZ_1120.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec24 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_3.IA,mean(FZ_1120.FZ)*ones(size(SA_vec)),tyre_coeffs);
MZ0_dfz_varGamma_vec25 = MF96_MZ_vec(zeros(size(SA_vec)),SA_vec,GAMMA_4.IA,mean(FZ_1120.FZ)*ones(size(SA_vec)),tyre_coeffs);

figure()
hold on
tiledlayout(2,3)
title('Self aligning moment - Raw data vs fitted curves - $\gamma = var$, $FZ = var$')
nexttile
hold on
plot(FZ_220.SA,FZ_220.MZ,'o')
plot(SA_vec,MZ0_dfz_varGamma_vec1)
plot(SA_vec,MZ0_dfz_varGamma_vec2)
plot(SA_vec,MZ0_dfz_varGamma_vec3)
plot(SA_vec,MZ0_dfz_varGamma_vec4)
plot(SA_vec,MZ0_dfz_varGamma_vec5)
title('FZ = 220')
legend('Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$')
xlabel('$\alpha$ [rad]')
ylabel('$M_{z} [N \cdot m]$')
nexttile
hold on
plot(FZ_440.SA,FZ_440.MZ,'o')
plot(SA_vec,MZ0_dfz_varGamma_vec6)
plot(SA_vec,MZ0_dfz_varGamma_vec7)
plot(SA_vec,MZ0_dfz_varGamma_vec8)
plot(SA_vec,MZ0_dfz_varGamma_vec9)
plot(SA_vec,MZ0_dfz_varGamma_vec10)
title('FZ = 440')
legend('Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$')
xlabel('$\alpha$ [rad]')
ylabel('$M_{z} [N \cdot m]$')
nexttile
hold on
plot(FZ_700.SA,FZ_700.MZ,'o')
plot(SA_vec,MZ0_dfz_varGamma_vec11)
plot(SA_vec,MZ0_dfz_varGamma_vec12)
plot(SA_vec,MZ0_dfz_varGamma_vec13)
plot(SA_vec,MZ0_dfz_varGamma_vec14)
plot(SA_vec,MZ0_dfz_varGamma_vec15)
title('FZ = 700')
legend('Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$')
xlabel('$\alpha$ [rad]')
ylabel('$M_{z} [N \cdot m]$')
nexttile
hold on
plot(FZ_900.SA,FZ_900.MZ,'o')
plot(SA_vec,MZ0_dfz_varGamma_vec16)
plot(SA_vec,MZ0_dfz_varGamma_vec17)
plot(SA_vec,MZ0_dfz_varGamma_vec18)
plot(SA_vec,MZ0_dfz_varGamma_vec19)
plot(SA_vec,MZ0_dfz_varGamma_vec20)
title('FZ = 900')
legend('Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$')
xlabel('$\alpha$ [rad]')
ylabel('$M_{z} [N \cdot m]$')
nexttile
hold on
plot(FZ_1120.SA,FZ_1120.MZ,'o')
plot(SA_vec,MZ0_dfz_varGamma_vec21)
plot(SA_vec,MZ0_dfz_varGamma_vec22)
plot(SA_vec,MZ0_dfz_varGamma_vec23)
plot(SA_vec,MZ0_dfz_varGamma_vec24)
plot(SA_vec,MZ0_dfz_varGamma_vec25)
title('FZ = 1120')
legend('Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$')
xlabel('$\alpha$ [rad]')
ylabel('$M_{z} [N \cdot m]$')


figure()
title('Self aligning moment - Raw data vs fitted curve - $\gamma = var$, $FZ = var$')
hold on
plot(ALPHA_vec,MZ_vec)
plot(SA_vec,MZ0_dfz_varGamma_vec1)
plot(SA_vec,MZ0_dfz_varGamma_vec2)
plot(SA_vec,MZ0_dfz_varGamma_vec3)
plot(SA_vec,MZ0_dfz_varGamma_vec4)
plot(SA_vec,MZ0_dfz_varGamma_vec5)
plot(SA_vec,MZ0_dfz_varGamma_vec6)
plot(SA_vec,MZ0_dfz_varGamma_vec7)
plot(SA_vec,MZ0_dfz_varGamma_vec8)
plot(SA_vec,MZ0_dfz_varGamma_vec9)
plot(SA_vec,MZ0_dfz_varGamma_vec10)
plot(SA_vec,MZ0_dfz_varGamma_vec11)
plot(SA_vec,MZ0_dfz_varGamma_vec12)
plot(SA_vec,MZ0_dfz_varGamma_vec13)
plot(SA_vec,MZ0_dfz_varGamma_vec14)
plot(SA_vec,MZ0_dfz_varGamma_vec15)
plot(SA_vec,MZ0_dfz_varGamma_vec16)
plot(SA_vec,MZ0_dfz_varGamma_vec17)
plot(SA_vec,MZ0_dfz_varGamma_vec18)
plot(SA_vec,MZ0_dfz_varGamma_vec19)
plot(SA_vec,MZ0_dfz_varGamma_vec20)
plot(SA_vec,MZ0_dfz_varGamma_vec21)
plot(SA_vec,MZ0_dfz_varGamma_vec22)
plot(SA_vec,MZ0_dfz_varGamma_vec23)
plot(SA_vec,MZ0_dfz_varGamma_vec24)
plot(SA_vec,MZ0_dfz_varGamma_vec25)
%% Longitudinal
%% Select tyre dataset
%dataset path
data_set_path = 'dataset/';
% dataset selection and loading

%data_set = 'Hoosier_B1464run23'; % pure lateral forces
data_set = 'Hoosier_B1464run30';  % braking/traction (pure log. force) + combined

% tyre geometric data:
% Hoosier	18.0x6.0-10
% 18 diameter in inches
% 6.0 section width in inches
% tread width in inches
diameter = 18*2.56; %
Fz0 = 220;   % [N] nominal load is given
R0  = diameter/2/100; % [m] get from nominal load R0 (m) *** TO BE CHANGED ***


fprintf('Loading dataset ...')
switch data_set
    case 'Hoosier_B1464run23'
  load ([data_set_path, 'Hoosier_B1464run23.mat']); % pure lateral
  cut_start = 27760;
  cut_end   = 54500;
    case 'Hoosier_B1464run30'
  load ([data_set_path, 'Hoosier_B1464run30.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 37643;
  otherwise 
  error('Not found dataset: `%s`\n', data_set) ;
  
end

% select dataset portion
smpl_range = cut_start:cut_end;

fprintf('completed!\n')
%% Plot raw data

figure
tiledlayout(6,1)

ax_list(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list,'x')

%% Select some specific data
% Cut crappy data and select only 12 psi data

vec_samples = 1:1:length(smpl_range);
tyre_data = table(); % create empty table
% store raw data in table
tyre_data.SL =  SL(smpl_range);
tyre_data.SA =  SA(smpl_range)*to_rad;
tyre_data.FZ = -FZ(smpl_range);  % 0.453592  lb/kg
tyre_data.FX =  FX(smpl_range);
tyre_data.FY =  FY(smpl_range);
tyre_data.MZ =  MZ(smpl_range);
tyre_data.IA =  IA(smpl_range)*to_rad;

% Extract points at constant inclination angle
GAMMA_tol = 0.05*to_rad;
idx.GAMMA_0 = 0.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 0.0*to_rad+GAMMA_tol;
idx.GAMMA_1 = 1.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 1.0*to_rad+GAMMA_tol;
idx.GAMMA_2 = 2.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 2.0*to_rad+GAMMA_tol;
idx.GAMMA_3 = 3.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 3.0*to_rad+GAMMA_tol;
idx.GAMMA_4 = 4.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 4.0*to_rad+GAMMA_tol;
idx.GAMMA_5 = 5.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 5.0*to_rad+GAMMA_tol;
GAMMA_0  = tyre_data( idx.GAMMA_0, : );
GAMMA_1  = tyre_data( idx.GAMMA_1, : );
GAMMA_2  = tyre_data( idx.GAMMA_2, : );
GAMMA_3  = tyre_data( idx.GAMMA_3, : );
GAMMA_4  = tyre_data( idx.GAMMA_4, : );
GAMMA_5  = tyre_data( idx.GAMMA_5, : );

% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

FZ_tol = 100;
idx.FZ_220  = 220-FZ_tol < tyre_data.FZ & tyre_data.FZ < 220+FZ_tol;
idx.FZ_440  = 440-FZ_tol < tyre_data.FZ & tyre_data.FZ < 440+FZ_tol;
idx.FZ_700  = 700-FZ_tol < tyre_data.FZ & tyre_data.FZ < 700+FZ_tol;
idx.FZ_900  = 900-FZ_tol < tyre_data.FZ & tyre_data.FZ < 900+FZ_tol;
idx.FZ_1120 = 1120-FZ_tol < tyre_data.FZ & tyre_data.FZ < 1120+FZ_tol;
FZ_220  = tyre_data( idx.FZ_220, : );
FZ_440  = tyre_data( idx.FZ_440, : );
FZ_700  = tyre_data( idx.FZ_700, : );
FZ_900  = tyre_data( idx.FZ_900, : );
FZ_1120 = tyre_data( idx.FZ_1120, : );

% The slip angle is varied continuously between -4 and +12° and then
% between -12° and +4° for the pure slip case

% The slip angle is varied step wise for longitudinal slip tests
% 0° , - 3° , -6 °
SA_tol = 0.5*to_rad;
idx.SA_0    =  0-SA_tol          < tyre_data.SA & tyre_data.SA < 0+SA_tol;
idx.SA_3neg = -(3*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -3*to_rad+SA_tol;
idx.SA_6neg = -(6*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -6*to_rad+SA_tol;
SA_0     = tyre_data( idx.SA_0, : );
SA_3neg  = tyre_data( idx.SA_3neg, : );
SA_6neg  = tyre_data( idx.SA_6neg, : );


figure()
tiledlayout(3,1)

ax_list(1) = nexttile;
plot(tyre_data.IA*to_deg)
hold on
plot(vec_samples(idx.GAMMA_0),GAMMA_0.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_1),GAMMA_1.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_2),GAMMA_2.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_3),GAMMA_3.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_4),GAMMA_4.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_5),GAMMA_5.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(2) = nexttile;
plot(tyre_data.FZ)
hold on
plot(vec_samples(idx.FZ_220),FZ_220.FZ,'.');
plot(vec_samples(idx.FZ_440),FZ_440.FZ,'.');
plot(vec_samples(idx.FZ_700),FZ_700.FZ,'.');
plot(vec_samples(idx.FZ_900),FZ_900.FZ,'.');
plot(vec_samples(idx.FZ_1120),FZ_1120.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')


ax_list(3) = nexttile;
plot(tyre_data.SA*to_deg)
hold on
plot(vec_samples(idx.SA_0),   SA_0.SA*to_deg,'.');
plot(vec_samples(idx.SA_3neg),SA_3neg.SA*to_deg,'.');
plot(vec_samples(idx.SA_6neg),SA_6neg.SA*to_deg,'.');
title('Slide slip')
xlabel('Samples [-]')
ylabel('[rad]')




%% Intersect tables to obtain specific sub-datasets

[TData0, ~] = intersect_table_data( SA_0, GAMMA_0, FZ_220 );

%% plot_selected_data

figure('Name','Selected-data')
plot_selected_data(TData0);


%% Fitting with Fz=Fz_nom= 220N and camber=0  alpha = 0 VX= 10
% ------------------
% long slip
tyre_coeffs.FZ0 = 220;
% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}
FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TData0.SL));
ones_vec  = ones(size(TData0.SL));

FX0_guess = MF96_FX0_vec(TData0.SL, zeros_vec, zeros_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% check guess 
figure('Name','FX0 guess')
plot(TData0.SL,TData0.FX,'.')
hold on
plot(TData0.SL,FX0_guess,'-')


% Plot raw data and initial guess
% figure()
% plot(TDataSub.KAPPA,TDataSub.FX,'o')
% hold on
% plot(TDataSub.KAPPA,FX0_guess,'x')

% Guess values for parameters to be optimised
%    [pCx1 pDx1 pEx1 pEx4  pHx1  pKx1  pVx1 
P0 = [  1,   2,   1,  0,   0,   1,   0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    [pCx1 pDx1 pEx1 pEx4  pHx1  pKx1  pVx1 
lb = [1,   0.1,   0,   0,  -10,    0,   -10];
ub = [2,    4,   1,   1,   10,   100,  10];


KAPPA_vec = TData0.SL;
FX_vec    = TData0.FX;

% check guess
SL_vec = -0.3:0.001:0.3;
FX0_fz_nom_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)) , zeros(size(SL_vec)), ...
                              FZ0.*ones(size(SL_vec)),tyre_coeffs);

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 

[P_fz_nom,fval,exitflag,iteration] = fmincon(@(P)resid_pure_Fx(P,FX_vec, KAPPA_vec,0,FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre data with new optimal values                             
tyre_coeffs.pCx1 = P_fz_nom(1) ; % 1
tyre_coeffs.pDx1 = P_fz_nom(2) ;  
tyre_coeffs.pEx1 = P_fz_nom(3) ;
tyre_coeffs.pEx4 = P_fz_nom(4) ;
tyre_coeffs.pHx1 = P_fz_nom(5) ; 
tyre_coeffs.pKx1 = P_fz_nom(6) ;
tyre_coeffs.pVx1 = P_fz_nom(7) ;

% R2 & RMSE
R2.FX0_FZnom = 1-resid_pure_Fx(P_fz_nom, FX_vec, KAPPA_vec,0,FZ0,tyre_coeffs);
RMSE.FX0_FZnom = sqrt((resid_pure_Fx(P_fz_nom, FX_vec, KAPPA_vec,0,FZ0,tyre_coeffs)*sum(FX_vec.^2))/length(FX_vec));

FX0_fz_nom_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)) , zeros(size(SL_vec)), ...
                              FZ0.*ones(size(SL_vec)),tyre_coeffs);
% Plots
figure('Name','Fx0(Fz0)')
plot(TData0.SL,TData0.FX,'o')
hold on
plot(SL_vec,FX0_fz_nom_vec,'-','LineWidth',2)
legend('Raw data', 'Fitted curve')
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')
title('Pure longitudinal force - Raw data vs fitted curve - FZ = 220 [N] $\gamma = 0$ [deg]')

%% Fit coefficient with variable load

% extract data with variable load
[TDataDFz, ~] = intersect_table_data( SA_0, GAMMA_0 );

zeros_vec = zeros(size(TDataDFz.SL));
ones_vec  = ones(size(TDataDFz.SL));
FX0_guess = MF96_FX0_vec(TDataDFz.SL,zeros_vec , zeros_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% check guess 
figure()
plot(TDataDFz.SL,TDataDFz.FX,'.')
hold on
plot(TDataDFz.SL,FX0_guess,'-')


% Guess values for parameters to be optimised
%    [pDx2 pEx2 pEx3 pHx2  pKx2  pKx3  pVx2] 
P0 = [  0,   0,   0,  0,   0,   0,   0]; 
lb = [];
ub = [];

KAPPA_vec = TDataDFz.SL;
FX_vec    = TDataDFz.FX;
FZ_vec    = TDataDFz.FZ;

% check guess
SL_vec = -0.3:0.001:0.3;
FX0_dfz_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)) , zeros(size(SL_vec)), ...
                           TDataDFz.FZ,tyre_coeffs);


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 

[P_dfz,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varFz(P,FX_vec, KAPPA_vec,0,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

disp(exitflag)
% Change tyre data with new optimal values                             
tyre_coeffs.pDx2 = P_dfz(1) ; % 1
tyre_coeffs.pEx2 = P_dfz(2) ;  
tyre_coeffs.pEx3 = P_dfz(3) ;
tyre_coeffs.pHx2 = P_dfz(4) ;
tyre_coeffs.pKx2 = P_dfz(5) ; 
tyre_coeffs.pKx3 = P_dfz(6) ;
tyre_coeffs.pVx2 = P_dfz(7) ;


res_FX0_dfz_vec = resid_pure_Fx_varFz(P_dfz,FX_vec,SL_vec,0 , FZ_vec,tyre_coeffs);


R2.FX0_FZvar = 1-resid_pure_Fx_varFz(P_dfz, FX_vec, KAPPA_vec,0,FZ_vec,tyre_coeffs);
RMSE.FX0_FZvar = sqrt((resid_pure_Fx_varFz(P_dfz, FX_vec, KAPPA_vec,0,FZ_vec,tyre_coeffs)*sum(FX_vec.^2))/length(FX_vec));

tmp_zeros = zeros(size(SL_vec));
tmp_ones = ones(size(SL_vec));


FX0_fz_var_vec1 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec2 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec3 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec4 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);


figure('Name','Fx0(Fz0)')
plot(TDataDFz.SL,TDataDFz.FX,'o')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
%plot(SL_vec,FX0_dfz_vec,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec1,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec2,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec3,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec4,'-','LineWidth',2)
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')
legend('Raw data','FZ = 220','FZ = 700', 'FZ = 900', 'FZ = 1120')
title('Pure longitudinal force - Raw data vs fitted curves - FZ = var [N] $\gamma = 0$ [deg]')


[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_220.FZ), tyre_coeffs);
Calfa_vec1_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_700.FZ), tyre_coeffs);
Calfa_vec2_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_900.FZ), tyre_coeffs);
Calfa_vec3_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_1120.FZ), tyre_coeffs);
Calfa_vec4_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);

Calfa_vec1 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec2 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec4 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

figure('Name','C_FK')

hold on

plot(SL_vec,Calfa_vec1,'-','LineWidth',2)
plot(SL_vec,Calfa_vec2,'-','LineWidth',2)
plot(SL_vec,Calfa_vec3,'-','LineWidth',2)
plot(SL_vec,Calfa_vec4,'-','LineWidth',2)
legend({'$Fz_{220}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})
xlabel('$\kappa$ [-]')
ylabel('$C_{F_\kappa}$ [N/-]')
title('Longitudinal slip stiffness')


%% Fit coefficient with variable camber

% extract data with variable load
[TDataGamma, ~] = intersect_table_data( SA_0, FZ_220 );

% Fit the coeffs { pDx3}

% Guess values for parameters to be optimised
P0 = [0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%lb = [0, 0,  0, 0,  0,  0,  0];
%ub = [2, 1e6,1, 1,1e1,1e2,1e2];
lb = [];
ub = [];


zeros_vec = zeros(size(TDataGamma.SL));
ones_vec  = ones(size(TDataGamma.SL));

KAPPA_vec = TDataGamma.SL;
GAMMA_vec = TDataGamma.IA; 
FX_vec    = TDataGamma.FX;
FZ_vec    = TDataGamma.FZ;

figure()
plot(KAPPA_vec,FX_vec);


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varGamma(P,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs.pDx3 = P_varGamma(1) ; % 1


R2.FX0_GAMMAvar = 1-resid_pure_Fx_varGamma(P_dfz, FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0,tyre_coeffs);
RMSE.FX0_GAMMAvar = sqrt((resid_pure_Fx_varGamma(P_dfz, FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0,tyre_coeffs)*sum(FX_vec.^2))/length(FX_vec));

FX0_varGamma_vec = MF96_FX0_vec(KAPPA_vec,zeros_vec , GAMMA_vec, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

figure('Name','Fx0 vs Gamma')
plot(KAPPA_vec,TDataGamma.FX,'o')
hold on
plot(KAPPA_vec,FX0_varGamma_vec,'-')
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')
title('Pure longitudinal force - Raw data vs fitted curve - FZ = 220 [N] $\gamma = var$ [deg]')
legend('Raw data','Fitted curve')
% Calculate the residuals with the optimal solution found above
res_Fx0_varGamma  = resid_pure_Fx_varGamma(P_varGamma,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_Fx0_varGamma);


[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs.FZ0, tyre_coeffs);
% 
fprintf('Bx      = %6.3f\n',Bx);
fprintf('Cx      = %6.3f\n',Cx);
fprintf('mux      = %6.3f\n',Dx/tyre_coeffs.FZ0);
fprintf('Ex      = %6.3f\n',Ex);
fprintf('SVx     = %6.3f\n',SVx);
fprintf('kappa_x = %6.3f\n',kappa__x);
fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_coeffs.FZ0);

%% Combined
%% Longitudinal fitting

[TData2,~] = intersect_table_data(FZ_220, GAMMA_0);
zeros_vec = zeros(size(TData2.SL));
ones_vec  = ones(size(TData2.SL));

KAPPA_vec = TData2.SL;
ALPHA_vec = TData2.SA; 
FX_vec    = TData2.FX;
FZ_vec    = TData2.FZ;

% Guess values for parameters to be optimised
%    [rBx1 rBx2 rCx1 rHx1] 
P0 = [ 8.3,  5,  0.9,   0];
lb = [   7,  0,  0.5,-100];
ub = [  20, 20,    3,   1];
[P_comb_x,fval,exitflag] = fmincon(@(P)resid_comb_Fx(P,FX_vec, KAPPA_vec, ALPHA_vec, zeros_vec,tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre_coeffs with new optimal values 
tyre_coeffs.rBx1 = P_comb_x(1);
tyre_coeffs.rBx2 = P_comb_x(2);
tyre_coeffs.rCx1 = P_comb_x(3);
tyre_coeffs.rHx1 = P_comb_x(4);

% R2 & RMSE
R2.FX_comb = 1-resid_comb_Fx(P_comb_x, FX_vec, KAPPA_vec, ALPHA_vec, zeros_vec,tyre_coeffs.FZ0, tyre_coeffs);
RMSE.FX_comb = sqrt((resid_comb_Fx(P_comb_x, FX_vec, KAPPA_vec, ALPHA_vec, zeros_vec,tyre_coeffs.FZ0, tyre_coeffs) ...
                            *sum(FX_vec.^2))/length(FX_vec));

% Plots
FX_comb = MF96_FX_vec(SL_vec,-3*to_rad*ones(length(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*tyre_coeffs.FZ0,tyre_coeffs);
FX_comb2 = MF96_FX_vec(SL_vec,-6*to_rad*ones(length(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*tyre_coeffs.FZ0,tyre_coeffs);

SA_vec = -0.3:0.001:0.3;
figure('Name','FX - Combined')
plot(KAPPA_vec,TData2.FX,'o')
hold on
plot(SL_vec,FX_comb,'-','LineWidth',2)
plot(SL_vec,FX_comb2,'-','LineWidth',2)
xlabel('$\kappa$ [-]')
ylabel('$F_{x}$ [N]')
legend('Raw data','$\alpha = -3 $ [deg]','$\alpha = -6 $ [deg]')
title('Combined longitudinal force - Raw data vs fitted curves')

SL_vec2 = -1:0.001:1;
[Gxa_veczero,~,~] = MF96_FXFYCOMB_vec(SL_vec2,zeros(length(SL_vec2)),...
                      zeros(length(SL_vec2)),ones(size(SL_vec2))*tyre_coeffs.FZ0,tyre_coeffs);
[Gxa_vecminus6,~,~] = MF96_FXFYCOMB_vec(SL_vec2,-6*to_rad*ones(length(SL_vec2)),...
                      zeros(length(SL_vec2)),ones(size(SL_vec2))*tyre_coeffs.FZ0,tyre_coeffs);
[Gxa_vecminus3,~,~] = MF96_FXFYCOMB_vec(SL_vec2,-3*to_rad*ones(length(SL_vec2)),...
                      zeros(length(SL_vec2)),ones(size(SL_vec2))*tyre_coeffs.FZ0,tyre_coeffs);
figure()
hold on
title('Combined longitudinal Gxa scaling factor')
plot(SL_vec2, Gxa_veczero(:,1),'LineWidth',2)
plot(SL_vec2, Gxa_vecminus3(:,1),'LineWidth',2)
plot(SL_vec2, Gxa_vecminus6(:,1),'LineWidth',2)
legend('$\alpha = 0$','$\alpha = -3 $ [deg]','$\alpha = -6$ [deg]')
xlabel('$\kappa$ [-]')
ylabel('$G_{xa}$ [-]')
%% Lateral fitting
[TData2,~] = intersect_table_data(FZ_900, GAMMA_0);

tyre_coeffs.FZ0 = 900;
FY_vec = TData2.FY;
KAPPA_vec = TData2.SL;
ALPHA_vec = TData2.SA; 
FX_vec    = TData2.FX;
FZ_vec    = TData2.FZ;
zeros_vec = zeros(size(TData2.SL));


% Guess values for parameters to be optimised
%    [rBy1 rBy2 rBy3 rCy1 rHy1 rVy1 rVy4 rVy5 rVy6] 
P0 = [ 4.9, 2.2,   0,   1,0.1,  0.1, 30,  0.5,  10];
lb = [-1000, -1000, -1000,-1000,-1000,-1000,-1000,-1000,-1000];
ub = [1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000];

options = optimset('MaxFunEvals',3.0e+10,'MaxIter',1.0e+5);
[P_comb_y,fval,exitflag] = fmincon(@(P)resid_comb_Fy(P,FY_vec, KAPPA_vec, ALPHA_vec, zeros_vec,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],[],[],[],options);

% Update tyre_coeffs with new optimal values 
tyre_coeffs.rBy1 = P_comb_y(1);
tyre_coeffs.rBy2 = P_comb_y(2);
tyre_coeffs.rBy3 = P_comb_y(3);
tyre_coeffs.rCy1 = P_comb_y(4);
tyre_coeffs.rHy1 = P_comb_y(5);
tyre_coeffs.rVy1 = P_comb_y(6);
tyre_coeffs.rVy4 = P_comb_y(7);
tyre_coeffs.rVy5 = P_comb_y(8);
tyre_coeffs.rVy6 = P_comb_y(9);
 
% R2 & RMSE
R2.FY_combFZnom = 1-resid_comb_Fy(P_comb_y,FY_vec, KAPPA_vec, ALPHA_vec, zeros_vec,FZ_vec, tyre_coeffs);
RMSE.FY_combFZnom = sqrt((resid_comb_Fy(P_comb_y,FY_vec, KAPPA_vec, ALPHA_vec, zeros_vec,FZ_vec, tyre_coeffs)*sum(FY_vec.^2))/length(FY_vec));

% Plots
FY_comb_0 = MF96_FY_vec(SL_vec,mean(SA_0.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*tyre_coeffs.FZ0,tyre_coeffs);
FY_comb_neg3 = MF96_FY_vec(SL_vec,mean(SA_3neg.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*tyre_coeffs.FZ0,tyre_coeffs);
FY_comb_neg6 = MF96_FY_vec(SL_vec,mean(SA_6neg.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*tyre_coeffs.FZ0,tyre_coeffs);

figure('Name','FY - Combined')
plot(KAPPA_vec,TData2.FY,'o')
hold on
plot(SA_vec,FY_comb_0,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg3,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg6,'-','LineWidth',2)
xlabel('$\kappa$ [-]')
ylabel('$F_{y}$ [N]')
title('Combined lateral force -Raw data vs fitted curves')
legend('Raw data','$\alpha = 0$ [deg]','$\alpha = -3$ [deg]','$\alpha = -6$ [deg]')

%% Lateral combined variable load
FY_vec = GAMMA_0.FY;
KAPPA_vec = GAMMA_0.SL;
ALPHA_vec = GAMMA_0.SA; 
FX_vec    = GAMMA_0.FX;
FZ_vec    = GAMMA_0.FZ;
zeros_vec = zeros(size(GAMMA_0.SL));


% Guess values for parameters to be optimised
%    [rVy2] 
P0 = [0];
lb = [];
ub = [];

options = optimset('MaxFunEvals',3.0e+10,'MaxIter',1.0e+5);
[P_comb_y_dfz,fval,exitflag] = fmincon(@(P)resid_comb_Fy_dfz(P,FY_vec, KAPPA_vec, ALPHA_vec, zeros_vec,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],[],[],[],options);

% Update tyre_coeffs with new optimal values 
tyre_coeffs.rVy2 = P_comb_y_dfz(1);


% R2 & RMSE
R2.FY_comb_dfz = 1-resid_comb_Fy_dfz(P_comb_y_dfz,FY_vec, KAPPA_vec, ALPHA_vec, zeros_vec,FZ_vec, tyre_coeffs);
RMSE.FY_comb_dfz = sqrt((resid_comb_Fy_dfz(P_comb_y_dfz,FY_vec, KAPPA_vec, ALPHA_vec, zeros_vec,FZ_vec, tyre_coeffs)*sum(FY_vec.^2))/length(FY_vec));

% Plots
FY_comb_0_1 = MF96_FY_vec(SL_vec,mean(SA_0.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*1120,tyre_coeffs);
FY_comb_neg3_1 = MF96_FY_vec(SL_vec,mean(SA_3neg.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*1120,tyre_coeffs);
FY_comb_neg6_1 = MF96_FY_vec(SL_vec,mean(SA_6neg.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*1120,tyre_coeffs);

FY_comb_0_2 = MF96_FY_vec(SL_vec,mean(SA_0.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*900,tyre_coeffs);
FY_comb_neg3_2 = MF96_FY_vec(SL_vec,mean(SA_3neg.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*900,tyre_coeffs);
FY_comb_neg6_2 = MF96_FY_vec(SL_vec,mean(SA_6neg.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*900,tyre_coeffs);

FY_comb_0_3 = MF96_FY_vec(SL_vec,mean(SA_0.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*700,tyre_coeffs);
FY_comb_neg3_3 = MF96_FY_vec(SL_vec,mean(SA_3neg.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*700,tyre_coeffs);
FY_comb_neg6_3 = MF96_FY_vec(SL_vec,mean(SA_6neg.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*700,tyre_coeffs);

FY_comb_0_4 = MF96_FY_vec(SL_vec,mean(SA_0.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*220,tyre_coeffs);
FY_comb_neg3_4 = MF96_FY_vec(SL_vec,mean(SA_3neg.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*220,tyre_coeffs);
FY_comb_neg6_4 = MF96_FY_vec(SL_vec,mean(SA_6neg.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*220,tyre_coeffs);

figure('Name','FY - Combined')
plot(KAPPA_vec,GAMMA_0.FY,'o')
hold on
plot(SA_vec,FY_comb_0_1,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg3_1,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg6_1,'-','LineWidth',2)
plot(SA_vec,FY_comb_0_2,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg3_2,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg6_2,'-','LineWidth',2)
plot(SA_vec,FY_comb_0_3,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg3_3,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg6_3,'-','LineWidth',2)
plot(SA_vec,FY_comb_0_4,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg3_4,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg6_4,'-','LineWidth',2)
xlabel('$\kappa$ [-]')
ylabel('$F_{y}$ [N]')
title('Combined lateral force -Raw data vs fitted curves')
legend('Raw data','$\alpha = 0$ [deg] $F_{Z} = 1120$ [N]','$\alpha = -3$ [deg] $F_{Z} = 1120$ [N]','$\alpha = -6$ [deg] $F_{Z} = 1120$ [N]' ...
    ,'$\alpha = 0$ [deg] $F_{Z} = 900$ [N]','$\alpha = -3$ [deg] $F_{Z} = 900$ [N]','$\alpha = -6$ [deg] $F_{Z} = 900$ [N]' ...
    ,'$\alpha = 0$ [deg] $F_{Z} = 700$ [N]','$\alpha = -3$ [deg] $F_{Z} = 700$ [N]','$\alpha = -6$ [deg] $F_{Z} = 700$ [N]' ...
    ,'$\alpha = 0$ [deg] $F_{Z} = 220$ [N]','$\alpha = -3$ [deg] $F_{Z} = 220$ [N]','$\alpha = -6$ [deg] $F_{Z} = 220$ [N]')
%% Lateral combined variable camber

FY_vec = FZ_900.FY;
KAPPA_vec = FZ_900.SL;
ALPHA_vec = FZ_900.SA; 
FX_vec    = FZ_900.FX;
FZ_vec    = FZ_900.FZ;
GAMMA_vec = FZ_900.IA;
zeros_vec = zeros(size(FZ_900.SL));


% Guess values for parameters to be optimised
%    [rVy3] 
P0 = [0];
lb = [];
ub = [];

options = optimset('MaxFunEvals',3.0e+10,'MaxIter',1.0e+5);
[P_comb_y_varGamma,fval,exitflag] = fmincon(@(P)resid_comb_Fy_varGamma(P,FY_vec, KAPPA_vec, ALPHA_vec, GAMMA_vec,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],[],[],[],options);

% Update tyre_coeffs with new optimal values 
tyre_coeffs.rVy3 = P_comb_y_varGamma(1);

% R2 & RMSE
R2.FY_comb_varGamma = 1-resid_comb_Fy_varGamma(P_comb_y,FY_vec, KAPPA_vec, ALPHA_vec, GAMMA_vec,FZ_vec, tyre_coeffs);
RMSE.FY_comb_varGamma = sqrt((resid_comb_Fy_varGamma(P_comb_y,FY_vec, KAPPA_vec, ALPHA_vec, GAMMA_vec,FZ_vec, tyre_coeffs)*sum(FY_vec.^2))/length(FY_vec));

%Plots
FY_comb_0_1 = MF96_FY_vec(SL_vec,mean(SA_0.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*900,tyre_coeffs);
FY_comb_neg3_1 = MF96_FY_vec(SL_vec,mean(SA_3neg.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*900,tyre_coeffs);
FY_comb_neg6_1 = MF96_FY_vec(SL_vec,mean(SA_6neg.SA)*ones(size(SL_vec)),zeros(length(SL_vec)),ones(length(SL_vec))*900,tyre_coeffs);

FY_comb_0_2 = MF96_FY_vec(SL_vec,mean(SA_0.SA)*ones(size(SL_vec)),2*to_rad*ones(length(SL_vec)),ones(length(SL_vec))*900,tyre_coeffs);
FY_comb_neg3_2 = MF96_FY_vec(SL_vec,mean(SA_3neg.SA)*ones(size(SL_vec)),2*to_rad*ones(length(SL_vec)),ones(length(SL_vec))*900,tyre_coeffs);
FY_comb_neg6_2 = MF96_FY_vec(SL_vec,mean(SA_6neg.SA)*ones(size(SL_vec)),2*to_rad*ones(length(SL_vec)),ones(length(SL_vec))*900,tyre_coeffs);

FY_comb_0_4 = MF96_FY_vec(SL_vec,mean(SA_0.SA)*ones(size(SL_vec)),4*to_rad*ones(length(SL_vec)),ones(length(SL_vec))*900,tyre_coeffs);
FY_comb_neg3_4 = MF96_FY_vec(SL_vec,mean(SA_3neg.SA)*ones(size(SL_vec)),4*to_rad*ones(length(SL_vec)),ones(length(SL_vec))*900,tyre_coeffs);
FY_comb_neg6_4 = MF96_FY_vec(SL_vec,mean(SA_6neg.SA)*ones(size(SL_vec)),4*to_rad*ones(length(SL_vec)),ones(length(SL_vec))*900,tyre_coeffs);

figure('Name','FY - Combined')
plot(KAPPA_vec,FZ_900.FY,'o')
hold on
plot(SA_vec,FY_comb_0_1,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg3_1,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg6_1,'-','LineWidth',2)
plot(SA_vec,FY_comb_0_2,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg3_2,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg6_2,'-','LineWidth',2)
plot(SA_vec,FY_comb_0_4,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg3_4,'-','LineWidth',2)
plot(SA_vec,FY_comb_neg6_4,'-','LineWidth',2)
xlabel('$\kappa$ [-]')
ylabel('$F_{y}$ [N]')
title('Combined lateral force -Raw data vs fitted curves')
legend('Raw data','$\alpha = 0$ [deg] $\gamma = 0$ [deg]','$\alpha = -3$ [deg] $\gamma = 0$ [deg]','$\alpha = -6$ [deg] $\gamma = 0$ [deg]' ...
    ,'$\alpha = 0$ [deg] $\gamma = 2$ [deg]','$\alpha = -3$ [deg] $\gamma = 2$ [deg]','$\alpha = -6$ [deg] $\gamma = 2$ [deg]' ...
    ,'$\alpha = 0$ [deg] $\gamma = 4$ [deg]','$\alpha = -3$ [deg] $\gamma = 4$ [deg]','$\alpha = -6$ [deg] $\gamma = 4$ [deg]')

[~,Gyk_vec,~] = MF96_FXFYCOMB_vec(SL_vec2,zeros(length(SL_vec2)),...
                      zeros(length(SL_vec2)),ones(size(SL_vec2))*tyre_coeffs.FZ0,tyre_coeffs);
[~,Gyk_vec_3neg,~] = MF96_FXFYCOMB_vec(SL_vec2,-6*to_rad*ones(length(SL_vec2)),...
                       zeros(length(SL_vec2)),ones(size(SL_vec2))*tyre_coeffs.FZ0,tyre_coeffs);
[~,Gyk_vec_6neg,~] = MF96_FXFYCOMB_vec(SL_vec2,-3*to_rad*ones(length(SL_vec2)),...
                       zeros(length(SL_vec2)),ones(size(SL_vec2))*tyre_coeffs.FZ0,tyre_coeffs);
figure()
hold on
plot(SL_vec2, Gyk_vec_3neg(:,1),'LineWidth',2)
plot(SL_vec2, Gyk_vec_6neg(:,1),'LineWidth',2)
title('Combined lateral $G_{yk}$ factor')
xlabel('$\kappa$ [-]')
ylabel('$G_{yk}$ [-]')
legend('$\alpha$ = -3 [deg]','$\alpha$ = -6 [deg]')




%% Save tyre data structure to mat file
%
save(['tyre_' data_set,'.mat'],'tyre_coeffs');


