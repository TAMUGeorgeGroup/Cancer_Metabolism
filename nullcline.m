% Function: nullcline
% Computes and plots the nullclines of a given mesh grid.
%
% Syntax:
%   out = nullcline(xmesh, ymesh, meshresult, varargin)
%
% Inputs:
%   xmesh      - Mesh grid for the x-axis values.
%   ymesh      - Mesh grid for the y-axis values.
%   meshresult - Matrix of computed values over the mesh grid.
%   varargin   - Additional arguments for contour plotting.
%
% Outputs:
%   out        - Contour matrix or handle to the contour plot.
%
% Description:
%   This function computes and plots the nullclines of a given mesh grid
%   using the contour or contourc functions. If four input arguments are
%   provided, the function uses contourc to compute the contour matrix.
%   Otherwise, it uses contour to plot the nullclines at the level [0 0].
%
% Example:
%   % Define mesh grids and computed values
%   [xmesh, ymesh] = meshgrid(0:0.1:10, 0:0.1:10);
%   meshresult = sin(xmesh) + cos(ymesh);
%
%   % Compute and plot nullclines
%   out = nullcline(xmesh, ymesh, meshresult, 'LineWidth', 2);
%
% See also:
%   contour, contourc

function out = nullcline(xmesh, ymesh, meshresult, varargin)
    % Check the number of input arguments
    if nargin == 4
        % Compute contour matrix using contourc
        out = contourc(xmesh(1,:), ymesh(:,1)', meshresult, varargin{1});
    else
        % Plot nullclines using contour at level [0 0]
        out = contour(xmesh, ymesh, meshresult, [0 0], varargin{:});
    end
end