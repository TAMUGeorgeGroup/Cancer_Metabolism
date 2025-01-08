% Function: findfixpoint_ggf
% Finds the fixed points of a given nullcline.
%
% Syntax:
%   [fixpoint, count] = findfixpoint_ggf(nullcline1, para1, para2, varargin)
%
% Inputs:
%   nullcline1 - Matrix containing the nullcline data.
%   para1      - Vector of parameters related to ROS and AMPK regulation.
%   para2      - Vector of parameters for the solve3d function.
%   varargin   - Additional arguments 
%
% Outputs:
%   fixpoint   - Matrix of fixed points (x, y values).
%   count      - Number of fixed points found.
%
% Description:
%   This function finds the fixed points of a given nullcline by solving a
%   system of nonlinear equations. It uses Hill functions to model the
%   interactions between key regulatory proteins and metabolic processes.
%   The function iterates over the nullcline data, solves the system of
%   equations, and identifies the fixed points where the sign of the
%   function changes.
%
%
% See also:
%   solve3d, Hillshift, Hillcomp

function [fixpoint, count] = findfixpoint_ggf(nullcline1, para1, para2, varargin)
    % Define Hill functions
    Hillshift = inline('lamda + (1 - lamda) * (1/(1 + (xx/xx0)^nxx))', 'xx', 'xx0', 'nxx', 'lamda');
    Hillcomp = inline('(gg0 + gg1*(xx1/xx10)^nxx1 + gg2*(xx2/xx20)^nxx2)/(1+(xx1/xx10)^nxx1+(xx2/xx20)^nxx2)', ...
                      'gg0', 'gg1', 'gg2', 'xx1', 'xx10', 'nxx1', 'xx2', 'xx20', 'nxx2');

    % Extract parameters for ROS and AMPK regulation
    gr = para1(1); 
    gama1 = para1(2);
    gamaf = para1(3);
    g0 = para1(4);
    g1 = para1(5);
    g2 = para1(6);
    Ar0_nox = para1(7);
    nar_nox = para1(8);
    hr0_nox = para1(9);
    nhr_nox = para1(10);
    gr_nox = para1(11);
    kr_nox = para1(12);
    ga = para1(13);   
    gh = para1(14);    
    Ah0 = para1(15); 
    Ra0 = para1(16); 
    Rh0 = para1(17);
    ha0 = para1(18); 
    G2h0 = para1(19); 
    nha = para1(20);  
    nah = para1(21);   
    nra = para1(22); 
    nrh = para1(23); 
    nG2h = para1(24); 
    kr = para1(25); 
    ka = para1(26);  
    kh = para1(27);
    Lamda_ra = para1(28); 
    Lamda_ha = para1(29); 
    Lamda_ah = para1(30);
    Lamda_rh = para1(31); 
    Lamda_G2h = para1(32);
    Aa0 = para1(33);
    naa = para1(34);
    Lamda_aa = para1(35);
    Ar0_n = para1(36);
    nar_n = para1(37);
    Lamda_ar_n = para1(38);
    gQ1 = para1(39);
    Q2R0 = para1(40);
    nQ2R = para1(41);
    Lamda_Q2R = para1(42);

    % Initialize count and fixpoint
    count = 0;
    fixpoint = [];
    xmax = varargin{1};

    % Extract nullcline data
    X = nullcline1(1, :);
    Y = nullcline1(2, :);
    
    % Iterate over nullcline data to find fixed points
    for index = 1:length(X)
        A = X(index);
        h = Y(index);
        
        % Solve the system of equations
        [x, G, C] = solve3d(A, h, para2);
        G1 = x(1);
        G2 = x(2);
        F = x(3);
        Q1 = x(4);
        Q3 = x(5);
        G3 = x(6);
        F2 = x(7);
        Q2 = x(8);
        
        % Calculate ATP production and ROS levels
        ATP = 2 * G2 + 29 * G1 + 106 * F + 24 * Q1 - 15 * Q3 - 13 * G3 - 7 * F2 - 2 * Q2;
        Rm = gr * (gama1 * G1 + gamaf * F + gQ1 * Q1) / (kr * (Hillshift(A, Ar0_n, nar_n, Lamda_ar_n) + 0.1 * Hillshift(Q2, Q2R0, nQ2R, Lamda_Q2R)));
        Rn = (gr_nox / (kr_nox * Hillshift(Q2, Q2R0, nQ2R, Lamda_Q2R))) * Hillcomp(g0, g1, g2, h, hr0_nox, nhr_nox, A, Ar0_nox, nar_nox);
        R = Rm + Rn;
        
        % Calculate the function value
        data_tmp = ga * Hillshift(R, Ra0, nra, Lamda_ra) * Hillshift(h, ha0, nha, Lamda_ha) * Hillshift(ATP, Aa0, naa, Lamda_aa) - ka * A;
        data(index) = data_tmp;
    end
    
    % Identify fixed points where the sign of the function changes
    for i = 1:(length(X) - 1)     
        if nullcline1(1, i) ~= 0 && nullcline1(1, i) ~= xmax && nullcline1(1, i + 1) ~= 0 && nullcline1(1, i + 1) ~= xmax
            Pre = data(i);
            Post = data(i + 1);
            if Pre * Post < 0 
                count = count + 1;
                fixpoint(count, 1) = X(i); % x value
                fixpoint(count, 2) = Y(i); % y value               
            end   
        end
    end
end