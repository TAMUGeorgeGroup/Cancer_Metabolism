% Function: root3d
% Computes the system of nonlinear equations representing metabolic pathways.
%
% Syntax:
%   Func = root3d(x, Parameters)
%
% Inputs:
%   x          - Vector of variables representing metabolic states.
%   Parameters - Vector of parameters required for the model.
%
% Outputs:
%   Func       - Vector of computed values for the system of equations.
%              - x = [G1, G2, F, Q1, Q3, G3, F2, Q2];  
%              - G1 (Glucose Oxidation), G2 (Glycolysis), F (Fatty Acid Oxidation), 
%                Q1 (Glutamine Oxidation), Q3 (Reductive Glutamine Metabolism), 
%                G3 (Reductive), Glucose Metabolism), F2 (Reductive Fatty Acid Metabolism), 
%                Q2 (Glutathione Synthesis), 
%
%
% Description:
%   This function defines a system of nonlinear equations that model the
%   dynamics of various metabolic pathways in cancer cells. The equations
%   are based on the interactions between key regulatory proteins and
%   metabolic processes. The function uses Hill functions to represent
%   these interactions.


function Func = root3d(x, Parameters)
    Hillshift = inline('lamda + (1 - lamda) * (1/(1 + (xx/xx0)^nxx))','xx', 'xx0', 'nxx','lamda');



    % Extract parameters from the input vector
    % A: AMPK level
    A = Parameters(1);
    
    % h: HIF-1 level
    h = Parameters(2);
    
    % g00: Basal glucose uptake rate
    g00 = Parameters(3);
    
    % c00: Basal acetyl-CoA utilization rate
    c00 = Parameters(4);
    
    % hg0: Threshold for HIF-1 regulation of glucose uptake
    hg0 = Parameters(5);
    
    % Ag0: Threshold for AMPK regulation of glucose uptake
    Ag0 = Parameters(6);
    
    % Ac0: Threshold for AMPK regulation of acetyl-CoA utilization
    Ac0 = Parameters(7);
    
    % nhg: Hill coefficient for HIF-1 regulation of glucose uptake
    nhg = Parameters(8);
    
    % nac: Hill coefficient for AMPK regulation of acetyl-CoA utilization
    nac = Parameters(9);
    
    % nag: Hill coefficient for AMPK regulation of glucose uptake
    nag = Parameters(10);
    
    % g1: Glucose oxidation rate
    g1 = Parameters(11);
    
    % g2: Glycolysis rate
    g2 = Parameters(12);
    
    % f0: Fatty acid oxidation rate
    f0 = Parameters(13);
    
    % ngg1: Hill coefficient for glucose oxidation regulation by glucose uptake
    ngg1 = Parameters(14);
    
    % ncg1: Hill coefficient for glucose oxidation regulation by acetyl-CoA utilization
    ncg1 = Parameters(15);
    
    % nhg2: Hill coefficient for HIF-1 regulation of glycolysis
    nhg2 = Parameters(16);
    
    % ngg2: Hill coefficient for glycolysis regulation by glucose uptake
    ngg2 = Parameters(17);
    
    % naf: Hill coefficient for AMPK regulation of fatty acid oxidation
    naf = Parameters(18);
    
    % ncf: Hill coefficient for fatty acid oxidation regulation by acetyl-CoA utilization
    ncf = Parameters(19);
    
    % hg2: Threshold for HIF-1 regulation of glycolysis
    hg2 = Parameters(20);
    
    % Af0: Threshold for AMPK regulation of fatty acid oxidation
    Af0 = Parameters(21);
    
    % lamda_hg: Baseline value for HIF-1 regulation of glucose uptake
    lamda_hg = Parameters(22);
    
    % lamda_ac: Baseline value for AMPK regulation of acetyl-CoA utilization
    lamda_ac = Parameters(23);
    
    % lamda_ag: Baseline value for AMPK regulation of glucose uptake
    lamda_ag = Parameters(24);
    
    % lamda_gg1: Baseline value for glucose oxidation regulation by glucose uptake
    lamda_gg1 = Parameters(25);
    
    % lamda_cg1: Baseline value for glucose oxidation regulation by acetyl-CoA utilization
    lamda_cg1 = Parameters(26);
    
    % lamda_hg2: Baseline value for HIF-1 regulation of glycolysis
    lamda_hg2 = Parameters(27);
    
    % lamda_gg2: Baseline value for glycolysis regulation by glucose uptake
    lamda_gg2 = Parameters(28);
    
    % lamda_af: Baseline value for AMPK regulation of fatty acid oxidation
    lamda_af = Parameters(29);
    
    % lamda_cf: Baseline value for fatty acid oxidation regulation by acetyl-CoA utilization
    lamda_cf = Parameters(30);
    
    % Fatty acid uptake parameters
    % gF00: Basal fatty acid uptake rate
    gF00 = Parameters(31);
    
    % AF00: Threshold for AMPK regulation of fatty acid uptake
    AF00 = Parameters(32);
    
    % nAF00: Hill coefficient for AMPK regulation of fatty acid uptake
    nAF00 = Parameters(33);
    
    % lamda_AF00: Baseline value for AMPK regulation of fatty acid uptake
    lamda_AF00 = Parameters(34);
    
    % nF00f1: Hill coefficient for fatty acid oxidation regulation by fatty acid uptake
    nF00f1 = Parameters(35);
    
    % lamda_F00f1: Baseline value for fatty acid oxidation regulation by fatty acid uptake
    lamda_F00f1 = Parameters(36);
    
    % HIF-1 repression of fatty acid oxidation
    % hf0: Threshold for HIF-1 repression of fatty acid oxidation
    hf0 = Parameters(37);
    
    % nhf: Hill coefficient for HIF-1 repression of fatty acid oxidation
    nhf = Parameters(38);
    
    % lamda_hf: Baseline value for HIF-1 repression of fatty acid oxidation
    lamda_hf = Parameters(39);
    
    % Glutamine oxidation parameters
    % q1: Glutamine oxidation rate
    q1 = Parameters(40);
    
    % hq0: Threshold for HIF-1 regulation of glutamine oxidation
    hq0 = Parameters(41);
    
    % nhq: Hill coefficient for HIF-1 regulation of glutamine oxidation
    nhq = Parameters(42);
    
    % lamda_hq: Baseline value for HIF-1 regulation of glutamine oxidation
    lamda_hq = Parameters(43);
    
    % Glutamine uptake parameters
    % nqq1: Hill coefficient for glutamine uptake regulation
    nqq1 = Parameters(44);
    
    % lamda_qq1: Baseline value for glutamine uptake regulation
    lamda_qq1 = Parameters(45);
    
    % q00: Basal glutamine uptake rate
    q00 = Parameters(46);
    
    % Reductive glutamine metabolism parameters
    % nqq3: Hill coefficient for reductive glutamine metabolism regulation
    nqq3 = Parameters(47);
    
    % lamda_qq3: Baseline value for reductive glutamine metabolism regulation
    lamda_qq3 = Parameters(48);
    
    % q3: Reductive glutamine metabolism rate
    q3 = Parameters(49);
    
    % HIF-1 and AMPK regulation of reductive glutamine metabolism
    % hq30: Threshold for HIF-1 regulation of reductive glutamine metabolism
    hq30 = Parameters(50);
    
    % nhq3: Hill coefficient for HIF-1 regulation of reductive glutamine metabolism
    nhq3 = Parameters(51);
    
    % lamda_hq3: Baseline value for HIF-1 regulation of reductive glutamine metabolism
    lamda_hq3 = Parameters(52);
    
    % Aq30: Threshold for AMPK regulation of reductive glutamine metabolism
    Aq30 = Parameters(53);
    
    % naq3: Hill coefficient for AMPK regulation of reductive glutamine metabolism
    naq3 = Parameters(54);
    
    % lamda_aq3: Baseline value for AMPK regulation of reductive glutamine metabolism
    lamda_aq3 = Parameters(55);
    
    % Reductive glucose metabolism parameters
    % ngg3: Hill coefficient for reductive glucose metabolism regulation
    ngg3 = Parameters(56);
    
    % lamda_gg3: Baseline value for reductive glucose metabolism regulation
    lamda_gg3 = Parameters(57);
    
    % g3: Reductive glucose metabolism rate
    g3 = Parameters(58);
    
    % AMPK repression of reductive glucose metabolism
    % Ag30: Threshold for AMPK repression of reductive glucose metabolism
    Ag30 = Parameters(59);
    
    % nag3: Hill coefficient for AMPK repression of reductive glucose metabolism
    nag3 = Parameters(60);
    
    % lamda_ag3: Baseline value for AMPK repression of reductive glucose metabolism
    lamda_ag3 = Parameters(61);
    
    % Reductive fatty acid metabolism parameters
    % gf2: Reductive fatty acid metabolism rate
    gf2 = Parameters(62);
    
    % nF00f2: Hill coefficient for reductive fatty acid metabolism regulation
    nF00f2 = Parameters(63);
    
    % lamda_F00f2: Baseline value for reductive fatty acid metabolism regulation
    lamda_F00f2 = Parameters(64);
    
    % Glutathione synthesis parameters
    % q2: Glutathione synthesis rate
    q2 = Parameters(65);
    
    % nqq2: Hill coefficient for glutathione synthesis regulation
    nqq2 = Parameters(66);
    
    % lamda_qq2: Baseline value for glutathione synthesis regulation
    lamda_qq2 = Parameters(67);
    
    % Part of G3 to F and part of Q3 to F
    % g3f: Fraction of G3 converted to F
    g3f = Parameters(68);
    
    % q3f: Fraction of Q3 converted to F
    q3f = Parameters(69);
    
    % MYC regulation parameters
    % Mg0: Threshold for MYC regulation of glucose uptake
    Mg0 = Parameters(70);
    
    % nMg: Hill coefficient for MYC regulation of glucose uptake
    nMg = Parameters(71);
    
    % lamda_Mg: Baseline value for MYC regulation of glucose uptake
    lamda_Mg = Parameters(72);
    
    % Mq0: Threshold for MYC regulation of glutamine uptake
    Mq0 = Parameters(73);
    
    % nMq: Hill coefficient for MYC regulation of glutamine uptake
    nMq = Parameters(74);
    
    % lamda_Mq: Baseline value for MYC regulation of glutamine uptake
    lamda_Mq = Parameters(75);
    
    % M: MYC level
    M = Parameters(76);
    
    % MYC regulation of reductive glucose metabolism
    % Mg30: Threshold for MYC regulation of reductive glucose metabolism
    Mg30 = Parameters(77);
    
    % nMg3: Hill coefficient for MYC regulation of reductive glucose metabolism
    nMg3 = Parameters(78);
    
    % lamda_Mg3: Baseline value for MYC regulation of reductive glucose metabolism
    lamda_Mg3 = Parameters(79);
    
    % MYC regulation of Q1
    % Mq10: Threshold for MYC regulation of Q1
    Mq10 = Parameters(80);
    
    % nMq1: Hill coefficient for MYC regulation of Q1
    nMq1 = Parameters(81);
    
    % lamda_Mq1: Baseline value for MYC regulation of Q1
    lamda_Mq1 = Parameters(82);
    
    % MYC regulation of Q2
    % Mq20: Threshold for MYC regulation of Q2
    Mq20 = Parameters(83);
    
    % nMq2: Hill coefficient for MYC regulation of Q2
    nMq2 = Parameters(84);
    
    % lamda_Mq2: Baseline value for MYC regulation of Q2
    lamda_Mq2 = Parameters(85);
    
    % MYC regulation of Q3
    % Mq30: Threshold for MYC regulation of Q3
    Mq30 = Parameters(86);
    
    % nMq3: Hill coefficient for MYC regulation of Q3
    nMq3 = Parameters(87);
    
    % lamda_Mq3: Baseline value for MYC regulation of Q3
    lamda_Mq3 = Parameters(88);
    
    % HIF-1 regulation of GSH synthesis
    % Aq20: Threshold for HIF-1 regulation of GSH synthesis
    Aq20 = Parameters(89);
    
    % nAq2: Hill coefficient for HIF-1 regulation of GSH synthesis
    nAq2 = Parameters(90);
    
    % lamda_Aq2: Baseline value for HIF-1 regulation of GSH synthesis
    lamda_Aq2 = Parameters(91);
    
    % HIF-1 activation of F2
    % hf20: Threshold for HIF-1 activation of F2
    hf20 = Parameters(92);
    
    % nhf2: Hill coefficient for HIF-1 activation of F2
    nhf2 = Parameters(93);
    
    % lamda_hf2: Baseline value for HIF-1 activation of F2
    lamda_hf2 = Parameters(94);
    
    
    %% Compute uptake rates using Hill functions
    % G0: Glucose uptake rate
    % C0: Acetyl-CoA utilization rate
    % Q0: Glutamine uptake rate
    % F00: Fatty acid uptake rate
    
    % Glucose uptake rate (G0)
    % G0 is influenced by HIF-1 (h), AMPK (A), and MYC (M)
    G0 = g00 * Hillshift(h, hg0, nhg, lamda_hg) + ...
         g00 * Hillshift(A, Ag0, nag, lamda_ag) + ...
         0 * g00 * Hillshift(M, Mg0, nMg, lamda_Mg);
    
    % Acetyl-CoA utilization rate (C0)
    % C0 is influenced by AMPK (A)
    C0 = c00 * Hillshift(A, Ac0, nac, lamda_ac);
    
    % Glutamine uptake rate (Q0)
    % Q0 is influenced by MYC (M)
    Q0 = q00 * Hillshift(M, Mq0, nMq, lamda_Mq);
    
    % Fatty acid uptake rate (F00)
    % F00 is influenced by AMPK (A)
    F00 = gF00 * Hillshift(A, AF00, nAF00, lamda_AF00);
    
    % System of equations for the metabolic pathways
    % Func: Vector of equations representing the metabolic states
    
    % Glucose oxidation (G1)
    % G1 is influenced by the total glucose uptake (G0) and acetyl-CoA utilization (C0)
    Func(1) = x(1) - g1 * Hillshift(x(1) + x(2) + x(6), G0, ngg1, lamda_gg1) * ...
                     Hillshift(2 * x(1) + 9 * x(3), C0, ncg1, lamda_cg1);
    
    % Glycolysis (G2)
    % G2 is influenced by HIF-1 (h) and the total glucose uptake (G0)
    Func(2) = x(2) - g2 * Hillshift(h, hg2, nhg2, lamda_hg2) * ...
                     Hillshift(x(1) + x(2) + x(6), G0, ngg2, lamda_gg2);
    
    % Fatty acid oxidation (F)
    % F is influenced by HIF-1 (h), fatty acid uptake (F00), AMPK (A), and acetyl-CoA utilization (C0)
    Func(3) = x(3) - f0 * Hillshift(h, hf0, nhf, lamda_hf) * ...
                     Hillshift(x(3) + x(7) - x(6) * g3f - x(5) * q3f, F00, nF00f1, lamda_F00f1) * ...
                     Hillshift(A, Af0, naf, lamda_af) * ...
                     Hillshift(2 * x(1) + 9 * x(3), C0, ncf, lamda_cf);
    
    % Glutamine oxidation (Q1)
    % Q1 is influenced by HIF-1 (h), the total glutamine uptake (Q0), and MYC (M)
    Func(4) = x(4) - q1 * Hillshift(h, hq0, nhq, lamda_hq) * ...
                     Hillshift(x(4) + x(5) + x(8), Q0, nqq1, lamda_qq1) * ...
                     Hillshift(M, Mq10, nMq1, lamda_Mq1);
    
    % Reductive glutamine metabolism (Q3)
    % Q3 is influenced by the total glutamine uptake (Q0), HIF-1 (h), AMPK (A), and MYC (M)
    Func(5) = x(5) - q3 * Hillshift(x(4) + x(5) + x(8), Q0, nqq3, lamda_qq3) * ...
                     Hillshift(h, hq30, nhq3, lamda_hq3) * ...
                     Hillshift(A, Aq30, naq3, lamda_aq3) * ...
                     Hillshift(M, Mq30, nMq3, lamda_Mq3);
    
    % Reductive glucose metabolism (G3)
    % G3 is influenced by the total glucose uptake (G0), AMPK (A), and MYC (M)
    Func(6) = x(6) - g3 * Hillshift(x(1) + x(2) + x(6), G0, ngg3, lamda_gg3) * ...
                     Hillshift(A, Ag30, nag3, lamda_ag3) * ...
                     Hillshift(M, Mg30, nMg3, lamda_Mg3);
    
    % Reductive fatty acid metabolism (F2)
    % F2 is influenced by the total fatty acid uptake (F00) and HIF-1 (h)
    Func(7) = x(7) - gf2 * Hillshift(x(3) + x(7) - x(6) * g3f - x(5) * q3f, F00, nF00f2, lamda_F00f2) * ...
                     Hillshift(h, hf20, nhf2, lamda_hf2);
    
    % Glutathione synthesis (Q2)
    % Q2 is influenced by the total glutamine uptake (Q0), AMPK (A), and MYC (M)
    Func(8) = x(8) - q2 * Hillshift(x(4) + x(5) + x(8), Q0, nqq2, lamda_qq2) * ...
                     Hillshift(A, Aq20, nAq2, lamda_Aq2) * ...
                     Hillshift(M, Mq20, nMq2, lamda_Mq2);
end
