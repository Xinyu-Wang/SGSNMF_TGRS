%  Default parameters of SGSNMF
%  Usage: para = default_SGSNMF()
%======================================================================
%          para.gama  -  regularization parameter; default{0.5} 
%           para.tol  -  tolerance for a relative stopping condition; default{0.05}
%       para.maxiter  -  Stop condition 1: limit of iterations; default{100}
%     para.timelimit  -  Stop condition 2: limit of time;
%       para.verbose  -  1- display current information in SGSNMF, and 0 - otherwise  
%    para.print_iter  -  iterations between print on screen
%======================================================================
% Copyright (C) 2018, Xinyu Wang (wangxinyu@whu.edu.cn)
%                     Yanfei Zhong (zhongyanfei@whu.edu.cn)
%                     Wuhan university
%                     All rights reserved.

% References:
% X. Wang, Y. Zhong, L. Zhang, and Y. Xu, ¡°Spatial Group Sparsity
% Regularized Nonnegative Matrix Factorization for Hyperspectral Unmixing,¡±
% IEEE Transactions on Geoscience and Remote Sensing, vol. 55, no. 11, pp.
% 6287-6304, NOV 2017, 2017.

% Last Modified:
% 27 Feb, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function para = default_SGSNMF()

para = {};
para.lambda = 0.4;
para.tol = 0.05;   
para.maxiter = 100; % stop condition 1
para.timelimit = 600; % stop condition 2
para.verbose = 1;
para.print_iter = 5;

return;