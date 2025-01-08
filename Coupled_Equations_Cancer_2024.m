%                  _      _ _                 ___                          
%  /\/\   ___   __| | ___| (_)_ __   __ _    / __\__ _ _ __   ___ ___ _ __ 
% /    \ / _ \ / _` |/ _ \ | | '_ \ / _` |  / /  / _` | '_ \ / __/ _ \ '__|
%/ /\/\ \ (_) | (_| |  __/ | | | | | (_| | / /__| (_| | | | | (_|  __/ |   
%\/    \/\___/ \__,_|\___|_|_|_| |_|\__, | \____/\__,_|_| |_|\___\___|_|   
%                                   |___/                                  

%{
The code contains the parameters representing cancer cells, it needs 
functions: 'solve3d' 'findfixpoint_ggf' 'root3d', 'nullcline'.
%}

clc; clear;

%% Functions
% Hillshift: Computes the Hill function with a shift.
% 
% Syntax:
%   result = Hillshift(xx, xx0, nxx, lamda)
%
% Inputs:
%   xx    - The variable of interest.
%   xx0   - The half-maximal effective concentration.
%   nxx   - The Hill coefficient.
%   lamda - The baseline value.
%


Hillshift = inline('lamda + (1 - lamda) * (1/(1 + (xx/xx0)^nxx))', 'xx', 'xx0', 'nxx', 'lamda');

% Hillcomp: Computes a composite Hill function.
% 
% Syntax:
%   result = Hillcomp(gg0, gg1, gg2, xx1, xx10, nxx1, xx2, xx20, nxx2)
%
% Inputs:
%   gg0   - Baseline value.
%   gg1   - Maximum effect of the first variable.
%   gg2   - Maximum effect of the second variable.
%   xx1   - The first variable of interest.
%   xx10  - The half-maximal effective concentration for the first variable.
%   nxx1  - The Hill coefficient for the first variable.
%   xx2   - The second variable of interest.
%   xx20  - The half-maximal effective concentration for the second variable.
%   nxx2  - The Hill coefficient for the second variable.
%
% Outputs:
%   result - The computed value of the composite Hill function.
%

Hillcomp = inline('(gg0 + gg1*(xx1/xx10)^nxx1 + gg2*(xx2/xx20)^nxx2)/(1+(xx1/xx10)^nxx1+(xx2/xx20)^nxx2)', ...
                  'gg0', 'gg1', 'gg2', 'xx1', 'xx10', 'nxx1', 'xx2', 'xx20', 'nxx2');






%% Parameters

% mtROS Level
% gr: Basal production rate of mtROS
% gama1: mtROS production in glucose oxidation
% gamaf: mtROS production in fatty acid oxidation
% kh: Degradation rate of HIF-1
% gQ1: mtROS production in glutamine oxidation
gr = 0.9;
gama1 = 50;
gamaf = gama1 * (9/2);
kh = 0.3;
gQ1 = 10;

% noxROS Level
% g0: Basal noxROS
% g1: Fold-change of noxROS activation by HIF-1
% g2: Fold-change of noxROS inhibition by AMPK
% Ar0_nox: Threshold of noxROS inhibition by AMPK
% nar_nox: Hill coefficient for noxROS inhibition by AMPK
% hr0_nox: Threshold of noxROS activation by HIF-1
% nhr_nox: Hill coefficient for noxROS activation by HIF-1
% gr_nox: Production rate of noxROS
% kr_nox: Degradation rate of noxROS
% ga: Production rate of AMPK
% gh: Production rate of HIF-1
% Ah0: Threshold of AMPK inhibition by HIF-1
% Ra0: Threshold of AMPK activation by ROS
% Rh0: Threshold of HIF-1 stabilization by ROS
% ha0: Threshold of HIF-1 activation by AMPK
% G2h0: Threshold of HIF-1 activation by glycolysis
g0 = 1.0;
g1 = 5;
g2 = 0.2;
Ar0_nox = 150;
nar_nox = 2;
hr0_nox = 250;
nhr_nox = 2;
gr_nox = 40;
kr_nox = 5.0;
ga = 40;   
gh = 15;    
Ah0 = 150; 
Ra0 = 250; 
Rh0 = 40;
ha0 = 150;  
G2h0 = 250;

% Hill Coefficients
% nha: Hill coefficient 
% nah: Hill coefficient 
% nra: Hill coefficient
% nrh: Hill coefficient 
% nG2h: Hill coefficient 
nha = 1;  
nah = 1;   
nra = 4; 
nrh = 4; 
nG2h = 4; 

% Degradation Rates
% kr: Degradation rate of mtROS
% ka: Degradation rate of AMPK
kr = 5.0; 
ka = 0.2;  

% Fold Changes
% Lamda_ra: Fold change of AMPK activation by ROS
% Lamda_ha: Fold change of AMPK inhibition by HIF-1
% Lamda_ah: Fold change of HIF-1 inhibition by AMPK
% Lamda_rh: Fold change of HIF-1 stabilization by ROS
% Lamda_G2h: Fold change of HIF-1 activation by glycolysis
Lamda_ra = 8;
Lamda_ha = 0.1; 
Lamda_ah = 0.1;
Lamda_rh = 0.2; 
Lamda_G2h = 0.1; 

% Parameters for AMPK Self-Inhibition
% Aa0: Threshold for AMPK inhibition by AMPK
% naa: Hill coefficient for AMPK inhibition
% Lamda_aa: Fold change for AMPK inhibition by ATP
Aa0 = 2000;
naa = 2;
Lamda_aa = 0.25;

% Parameters for AMPK Negative Regulation of mtROS Pathway
% Ar0_n: Threshold of mtROS inhibition by AMPK
% nar_n: Hill coefficient 
% Lamda_ar_n: Fold change for AMPK regulation of mtROS pathway
Ar0_n = 350;
nar_n = 2;
Lamda_ar_n = 2;  

% GSH on ROS
% Q2R0: Threshold for GSH regulation of ROS
% nQ2R: Hill coefficient for GSH regulation of ROS
% Lamda_Q2R: Fold change of mtROS clearance by AMPK
Q2R0 = 10;
nQ2R = 2;
Lamda_Q2R = 1.2;

% MYC Parameters
% M: MYC level
% Mh0: Threshold of MYC inhibting HIF-1 degradation
% nMh: Hill coefficient 
% Lamda_Mh: Fold change of MYC inhibting HIF-1 degradation
M = 1200;
Mh0 = 150;
nMh = 2;
Lamda_Mh = 0.8;

% MYC Promotes Glucose Uptake
% Mg0: Threshold of glucose uptake regulated by MYC
% nMg: Hill coefficient 
% lamda_Mg: Fold change of glucose uptake regulated by MYC
Mg0 = 150;
nMg = 2;
lamda_Mg = 2;

% MYC Promotes Glutamine Uptake
% Threshold of glutamine uptake regulated by MYC
% nMq: Hill coefficient 
% lamda_Mq: Fold-change of glutamine uptake regulated by MYC
Mq0 = 150; 
nMq = 2;
lamda_Mq = 2; 

% MYC Reductive Glucose Metabolism
% Mg30: Threshold of reductive glucose metabolism regulation by MYC
% nMg3: Hill coefficient 
% lamda_Mg3: Fold-change of reductive glucose metabolism regulation by MYC
Mg30 = 150; 
nMg3 = 2;
lamda_Mg3 = 2;

% MYC Regulation of Q1
% Mq10: Threshold of glutamine oxidation regulation by MYC
% nMq1: Hill coefficient for MYC regulation of Q1
% lamda_Mq1: Fold change of glutamine oxidation regulation by MYC
Mq10 = 150;
nMq1 = 2;
lamda_Mq1 = 2; 

% MYC Regulation of Q2
% Mq20: Threshold of glutathione synthesis regulation by MYC
% nMq2: Hill coefficient for MYC regulation of Q2
% lamda_Mq2: Fold change of glutathione synthesis regulation by MYC
Mq20 = 150;
nMq2 = 2;
lamda_Mq2 = 0.9; 

% MYC Regulation of Q3
% Mq30: Threshold of reductive glutamine metabolism regulation by MYC
% nMq3: Hill coefficient for MYC regulation of Q3
% lamda_Mq3: Fold change of reductive glutamine metabolism regulation by MYC
Mq30 = 150;
nMq3 = 2;
lamda_Mq3 = 2; 

% Combine Parameters into a Vector
para1 = [gr, gama1, gamaf, g0, g1, g2, Ar0_nox, nar_nox, hr0_nox, nhr_nox, gr_nox, kr_nox, ga, gh, ...
         Ah0, Ra0, Rh0, ha0, G2h0, nha, nah, nra, nrh, nG2h, kr, ka, kh, Lamda_ra, Lamda_ha, ...
         Lamda_ah, Lamda_rh, Lamda_G2h, Aa0, naa, Lamda_aa, Ar0_n, nar_n, Lamda_ar_n, gQ1, ...
         Q2R0, nQ2R, Lamda_Q2R];

%% Parameters for the Metabolic Pathways

% Glucose Uptake and Utilization
% g00: Basal glucose uptake rate
% c00: Basal acetyl-CoA utilization rate
% hg0: Threshold for HIF-1 regulation of glucose uptake
% Ag0: Threshold for AMPK regulation of glucose uptake
% Ac0: Threshold for AMPK regulation of acetyl-CoA utilization
g00 = 20;
c00 = 15; 
hg0 = 150;
Ag0 = 200;
Ac0 = 250;

% Basal Rates for Metabolic Pathways
% g1_base: Basal glucose oxidation rate
% g2_base: Basal glycolysis rate
% f0: Basal fatty acid oxidation rate
g1_base = 100;
g2_base = 150;
f0 = 3;

% Thresholds for HIF-1 and AMPK Regulation
% hg2: Threshold of glycolysis upregulation by HIF-1
% Af0: Threshold for AMPK regulation of fatty acid oxidation
hg2 = 200;
Af0 = 200;

% Hill Coefficients for Metabolic Pathways
% nhg: Hill coefficient for HIF-1 regulation of glucose uptake
% nac: Hill coefficient for AMPK regulation of acetyl-CoA utilization
% nag: Hill coefficient for AMPK regulation of glucose uptake
% ncg1: Hill coefficient for glucose oxidation regulation by acetyl-CoA utilization
% nhg2: Hill coefficient for HIF-1 regulation of glycolysis
% ngg1: Hill coefficient for glucose oxidation regulation by glucose uptake
% ngg2: Hill coefficient for glycolysis regulation by glucose uptake
% naf: Hill coefficient for AMPK regulation of fatty acid oxidation
% ncf: Hill coefficient for fatty acid oxidation regulation by acetyl-CoA utilization
nhg = 4;
nac = 4;
nag = 2;
ncg1 = 4;
nhg2 = 4;
ngg1 = 2;
ngg2 = 2;
naf = 4;
ncf = 2;

% Baseline Values for Hill Functions
% lamda_hg: Fold-change of glucose uptake regulated by HIF-1
% lamda_ag: Fold-change of glucose uptake regulated by AMPK
% lamda_ac: Fold-change of acetyl-CoA utilization regulated by AMPK
% lamda_cg1: Restriction of glucose oxidation by acetyl-CoA utilization rate
% lamda_gg1: Restriction of glucose oxidation rate by glucose uptake rate
% lamda_gg2: Restriction of glycolysis rate by glucose uptake rate
% lamda_cf: Restriction of FAO rate by acetyl-CoA utilization rate
% lamda_af: Fold-change of FAO upregulation by AMPK
% lamda_hg2: Fold-change of glycolysis upregulation by HIF-1
lamda_hg = 6;
lamda_ag = 4;
lamda_ac = 8;
lamda_cg1 = 0.1;
lamda_gg1 = 0.1;
lamda_gg2 = 0.1;
lamda_cf = 0.1;
lamda_af = 6;
lamda_hg2 = 8;

% Fatty Acid Uptake Parameters
% gF00: Basal fatty acid uptake rate
% AF00: Threshold for AMPK regulation of fatty acid uptake
% nAF00: Hill coefficient for AMPK regulation of fatty acid uptake
% lamda_AF00: Fold-change of fatty acid uptake regulation by AMPK
% nF00f1: Hill coefficient for fatty acid oxidation regulation by fatty acid uptake
% lamda_F00f1: Restriction of FAO rate by fatty acid uptake rate
gF00 = 5;
AF00 = 200;
nAF00 = 2;
lamda_AF00 = 4;
nF00f1 = 2;
lamda_F00f1 = 0.8;

% HIF-1 Repression of Fatty Acid Oxidation (FAO)
% hf0: Threshold for HIF-1 regulation of FAO
% nhf: Hill coefficient for HIF-1 repression of FAO
% lamda_hf: Fold-change of FAO regulation by HIF-1
hf0 = 200;
nhf = 4;
lamda_hf = 0.8;

% Glutamine Metabolism Parameters
% q1: Glutamine oxidation rate
% hq0: Threshold for HIF-1 regulation of glutamine oxidation
% nhq: Hill coefficient for HIF-1 regulation of glutamine oxidation
% lamda_hq: Fold change of glutamine oxidation regulation by HIF-1
% nqq1: Hill coefficient for glutamine uptake regulation
% lamda_qq1: Restriction of glutamine oxidation rate by glutamine uptake rate
% q00: Basal glutamine uptake rate
q1 = 30;
hq0 = 350;
nhq = 4;
lamda_hq = 0.2;
nqq1 = 2;
lamda_qq1 = 0.5;
q00 = 10;

% Reductive Glutamine Metabolism Parameters (Q3)
% nqq3: Hill coefficient for reductive glutamine metabolism regulation
% lamda_qq3: Baseline value for reductive glutamine metabolism regulation
% q3: Basal reductive glutamine metabolic rate
nqq3 = 2;
lamda_qq3 = 0.5;
q3 = 20;

% HIF-1 and AMPK Regulation of Reductive Glutamine Metabolism
% hq30: Threshold for HIF-1 regulation of reductive glutamine metabolism
% nhq3: Hill coefficient for HIF-1 regulation of reductive glutamine metabolism
% lamda_hq3: Fold change of reductive glutamine metabolism regulation by HIF-1
% Aq30: Threshold of reductive glutamine metabolism regulation by AMPK
% naq3: Hill coefficient for AMPK regulation of reductive glutamine metabolism
% lamda_aq3: Fold change of reductive glutamine metabolism regulation by AMPK
hq30 = 200;
nhq3 = 4;
lamda_hq3 = 2;
Aq30 = 200;
naq3 = 2;
lamda_aq3 = 0.5;

% Reductive Glucose Metabolism Parameters
% ngg3: Hill coefficient for reductive glucose metabolism regulation
% lamda_gg3: Restriction of reductive glucose metabolic rate by glucose uptake rate
% g3: Reductive glucose metabolism rate
ngg3 = 2;
lamda_gg3 = 0.1;
g3 = 60;

% AMPK Repression of Reductive Glucose Metabolism
% Ag30: Threshold of reductive glucose metabolism regulation by AMPK
% nag3: Hill coefficient for AMPK repression of reductive glucose metabolism
% lamda_ag3: Fold-change of reductive glucose metabolism regulation by AMPK
Ag30 = 200;
nag3 = 2;
lamda_ag3 = 0.25;

% Reductive Fatty Acid Metabolism Parameters
% gf2: Reductive fatty acid metabolism rate
% nF00f2: Hill coefficient for reductive fatty acid metabolism regulation
% lamda_F00f2: Restriction of FAO rate by fatty acid uptake rate
gf2 = 75;
nF00f2 = 2;
lamda_F00f2 = 0.1;

% Glutathione (GSH) Production Parameters (Q2)
% q2: GSH production rate
% nqq2: Hill coefficient for GSH production regulation
% lamda_qq2: Restriction of glutathione synthesis rate by glutamine uptake rate
q2 = 30;
nqq2 = 2;
lamda_qq2 = 0.5;

% Conversion Factors for G3 and Q3 to F
% g3f: Lipogenesis rate of glucose
% q3f: Lipogenesis rate of glutamine
g3f = 0.1;
q3f = 0.1;

% AMPK Regulation of GSH
% Aq20: Threshold for AMPK regulation of GSH
% nAq2: Hill coefficient for AMPK regulation of GSH
% lamda_Aq2: Baseline value for AMPK regulation of GSH
Aq20 = 250;
nAq2 = 4;
lamda_Aq2 = 5;

% HIF-1 Activation of F2
% hf20: Threshold for HIF-1 activation of F2
% nhf2: Hill coefficient for HIF-1 activation of F2
% lamda_hf2: Baseline value for HIF-1 activation of F2
hf20 = 200;
nhf2 = 2;
lamda_hf2 = 4;

% Combine Parameters into a Vector
para = [g00, c00, hg0, Ag0, Ac0, nhg, nac, nag, g1_base, g2_base, f0, ngg1, ncg1,...
        nhg2, ngg2, naf, ncf, hg2, Af0, lamda_hg, lamda_ac, lamda_ag, ...
        lamda_gg1, lamda_cg1, lamda_hg2, lamda_gg2, lamda_af, lamda_cf, ...
        gF00, AF00, nAF00, lamda_AF00, nF00f1, lamda_F00f1, hf0, nhf, lamda_hf, q1, ...
        hq0, nhq, lamda_hq, nqq1, lamda_qq1, q00, nqq3, lamda_qq3, q3, hq30, nhq3, ...
        lamda_hq3, Aq30, naq3, lamda_aq3, ngg3, lamda_gg3, g3, Ag30, nag3, lamda_ag3, ... 
        gf2, nF00f2, lamda_F00f2, q2, nqq2, lamda_qq2, g3f, q3f, Mg0, nMg, lamda_Mg, Mq0, ...
        nMq, lamda_Mq, M, Mg30, nMg3, lamda_Mg3, Mq10, nMq1, lamda_Mq1, Mq20, nMq2, ...
        lamda_Mq2, Mq30, nMq3, lamda_Mq3, Aq20, nAq2, lamda_Aq2, hf20, nhf2, lamda_hf2];



%% Nullclines
% Define the range and resolution for AMPK and HIF-1 levels
r = 750;
N = 50;
xmax = r;
ymax = r;

dx = xmax / N;
dy = ymax / N;

A_vec = 0:dx:xmax;
h_vec = 0:dy:ymax;

% Timer
tic

% Compute nullclines
for i = 1:length(A_vec)
    % Print progress
    fprintf('Iteration %d of %d\n', i, length(A_vec));

    A = A_vec(i);
    for j = 1:length(h_vec)
        h = h_vec(j);
        % Solve the system of equations for given A and h
        [x, G, C] = solve3d(A, h, para);
        G1 = x(1);
        G2 = x(2);
        F = x(3);
        Q1 = x(4);
        Q3 = x(5);
        G3 = x(6);
        F2 = x(7);
        Q2 = x(8);

        % Calculate ATP production
        ATP = 2 * G2 + 29 * G1 + 106 * F + 24 * Q1 - 15 * Q3 - 13 * G3 - 7 * F2 - 2 * Q2;
        G = G1 + G2 + G3;
        C = 2 * G1 + 9 * F;

        % Effect of GSH on ROS
        Rm = gr * (gama1 * G1 + gamaf * F + gQ1 * Q1) / (kr * (Hillshift(A, Ar0_n, nar_n, Lamda_ar_n) + 0.1 * Hillshift(Q2, Q2R0, nQ2R, Lamda_Q2R)));
        Rn = (gr_nox / (kr_nox * Hillshift(Q2, Q2R0, nQ2R, Lamda_Q2R))) * Hillcomp(g0, g1, g2, h, hr0_nox, nhr_nox, A, Ar0_nox, nar_nox);
        R = Rm + Rn;

        % Calculate the equations for AMPK and HIF-1
        eqsA = ga * Hillshift(R, Ra0, nra, Lamda_ra) * Hillshift(h, ha0, nha, Lamda_ha) * Hillshift(ATP, Aa0, naa, Lamda_aa) - ka * A;
        eqsh = gh * Hillshift(A, Ah0, nah, Lamda_ah) - kh * h * Hillshift(R, Rh0, nrh, Lamda_rh) * Hillshift(G2, G2h0, nG2h, Lamda_G2h) * Hillshift(M, Mh0, nMh, Lamda_Mh);
        
        result1(j, i) = eqsA;
        result2(j, i) = eqsh;
    end
end

%% Plot the nullclines and find the steady states
figure1 = figure('Color', [1 1 1]);
axes1 = axes('Parent', figure1, ...
    'Position', [0.20323275862069 0.17109634551495 0.570474137931034 0.778820598006645], ...
    'PlotBoxAspectRatio', [1 1 1], ...
    'LineWidth', 2, ...
    'FontSize', 30);
hold(axes1, 'all');
box on;
xlim([0, r]);
ylim([0, r]);

% Plot nullclines
Ag_nullcline = nullcline(A_vec, h_vec, result1, 'LineWidth', 4, 'LineColor', [0 0 1], 'LevelList', 0);
hg_nullcline = nullcline(A_vec, h_vec, result2, 'LineWidth', 4, 'LineColor', [1 0 0], 'LevelList', 0);

% Find fixed points
[fixpoint, count] = findfixpoint_ggf(hg_nullcline, para1, para, xmax);

N_init = size(fixpoint, 1);
atp_production = [];
stable = [];
for i = 1:N_init
    A = fixpoint(i, 1);
    h = fixpoint(i, 2);
    [x, G, C] = solve3d(A, h, para);
    G1 = x(1);
    G2 = x(2);
    F = x(3);
    Q1 = x(4);
    Q3 = x(5);
    G3 = x(6);
    F2 = x(7);
    Q2 = x(8);

    % Calculate ROS levels
    Rm = gr * (gama1 * G1 + gamaf * F + gQ1 * Q1) / (kr * (Hillshift(A, Ar0_n, nar_n, Lamda_ar_n) + 0.1 * Hillshift(Q2, Q2R0, nQ2R, Lamda_Q2R)));
    Rn = (gr_nox / (kr_nox * Hillshift(Q2, Q2R0, nQ2R, Lamda_Q2R))) * Hillcomp(g0, g1, g2, h, hr0_nox, nhr_nox, A, Ar0_nox, nar_nox);
    
    % Store ATP production and stable states
    atp_production = [atp_production; [29 * G1, 2 * G2, 106 * F, 24 * Q1, -15 * Q3, -13 * G3, -7 * F2, -2 * Q2]];
    stable = [stable; [A, h, Rm, Rn, G1, G2, G3, F, F2, Q1, Q2, Q3]];
end

%% Plot the data
hold on;
if size(stable, 1) ~= 0
    X1 = stable(1:2:5, 1);
    Y1 = stable(1:2:5, 2);
    h3 = plot(X1, Y1, 'or', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'g');
    X2 = stable(2:2:4, 1);
    Y2 = stable(2:2:4, 2);
    h4 = plot(X2, Y2, 'or', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'w');
end

xlabel('AMPK');
ylabel('HIF1');
legend('pAMPK', 'HIF1', 'Stable', 'Unstable');
hold off;

% Define names for rows and variables (columns)
row_names = {'O', 'W/O', 'W'};
column_names = {'G1', 'G2', 'F', 'Q1', 'Q3', 'G3', 'F2', 'Q2'};

% Convert array to table
atp_production = array2table(atp_production);

% Assign names to rows and variables
atp_production.Properties.VariableNames = column_names;

% Extract data from table
data = table2array(atp_production);

% Create bar plot
figure;
bar(data(1:2:end, :));
xlabel('Columns');
ylabel_handle = ylabel('ATP production rate (µM/s)');
set(ylabel_handle, 'FontSize', 20); % Increase font size for ylabel

% Add legend
legend(column_names);

% Increase font size of axis labels
set(gca, 'FontSize', 14); % Adjust the font size as needed

toc