% Matlab Demo of SGSNMF
%
% Copyright (C) 2018, Xinyu Wang (wangxinyu@whu.edu.cn)
%                     Yanfei Zhong (zhongyanfei@whu.edu.cn)
%                     Wuhan university
%                     All rights reserved.
%
% References:
% X. Wang, Y. Zhong, L. Zhang, and Y. Xu, "Spatial Group Sparsity
% Regularized Nonnegative Matrix Factorization for Hyperspectral Unmixing",
% IEEE Transactions on Geoscience and Remote Sensing, vol. 55, no. 11, pp.
% 6287-6304, 2017.
%
% Notes:
%
%    If a higher accuracy is needed, some parameters, such as Sw(the
%    averange width of superpixels), Ws(Trade-off coefficient between
%    spatial and spectral distances) and para.lambda(Trade-off coefficient
%    between the reconstruction and the group sparsity regularizer) can be
%    adjusted by the user referring to the given range and the HSI.
%
% Last Modified:
% 27 Dem, 2018
% ********************************************************************** %

close all
clear all
clc
addpath(genpath('code'));
addpath(genpath('../data'));
 
% ********************************************************************** %
% load the groundtruth endmember and abudaces matrix
load 'F1_A_9'; % A L*M endmember matrix
load 'F1_S_9'; % S row*col*M abundance maps

[row, col, M] = size(S);
N = row*col;
L = size(A,1); 
S = reshape(S,N,M)'; % colomnwise abundance matrix

% build 3D noiseless HSI
Y = reshape((A*S)',row,col,L);  % For real data, please normalize the radiance to a range of 0 - 1.0.

% add white Gaussian noise
SNR = 20;
noise_type = 'additive'; eta = 0;
[X, n, Cn] = addNoise (Y,noise_type,SNR, eta, 1);

% nonnegative obervation
X = max(X,eps); % 2D HSI L*N
Y = reshape(X',row,col,L);  % 3D HSI row*col*L

%% ********************* Pre-processing ********************** %%
% ## 1.SLIC image segmentation ****************************************** %
% Sw - the averange width of superpixels; {3-11} 
% P  - the number of superpixels;
% Ws - Trade-off coefficient between spatial and spectral distances;{0.5}
Sw = 6; P = round(row*col/Sw^2); Ws = 0.5;
seg = slic_HSI(Y, P, Ws);

% ## 2.Initialize A and S *********************************************** %
% use region-based vca to initiate A_init
[A_VCA,~] = hyperVca(X,M);
[A_init,~] = hyperVca(seg.X_c,M);
% use FCLS to initiate S_init
S_init = fcls(A_init,X);

% ## plot the initial results ******************************************* %
% show the initial accuracy(SAD and RMSE)
fprintf('initial SAD estimation:\n');
Sam = sam(A, A_init); 
fprintf('initial RMSE estimation:\n');
rmse(S, S_init, Sam(1,:), Sam(2,:));

figure;
p_row = ceil((M+1)./3); 
for i = 1 : M        
    subplot_tight(3, p_row, i,[.01 .01]);imagesc(reshape(S_init(i,:)',row, col),[0,1]);axis image;axis off;
end
subplot_tight(3, p_row, M+1,[.01 .01]);imagesc(reshape(sqrt(sum((X - A_init*S_init).^2))',row, col),[0,1]);axis image;axis off;
drawnow;
saveas(gcf,'../results/initial_abun.tif');

%% *************** SGSNMF *************** %%
% ## 3.main function of SGSNMF ****************************************** %
% show the parameters and initial matrices 
para = default_SGSNMF;
para.Y =Y; 
para.X = X; 
para.M = M; 
para.W =A_init; 
para.H =S_init;
para.lambda = 0.3;
fprintf('\n parameter settings: \n');
disp(para); disp(seg);

% W - estimated endmember martix;
% H - estimated abundance matrix;
[W,H] = sgsnmf(para,seg); 

% ## plot the finial results ****************************************** %
% show SAD and RMSE results of SGSNMF
fprintf('final_SAD estimated: \n');
Sam = sam(A, W); 
fprintf('final_RMSE estimated: \n');
rmse(S, H, Sam(1,:), Sam(2,:));

%  plot estimated abundances maps
figure;
for i = 1 : para.M        
    subplot_tight(3, p_row, i,[.01 .01]);
    imagesc(reshape(H(i,:)',row, col),[0,1]);axis image;axis off;
end
subplot_tight(3, p_row, para.M+1,[.01 .01]);
imagesc(reshape(sqrt(sum((X-W*H).^2))',row, col),[0,1]);axis image;axis off;
drawnow;
saveas(gcf,'../results/estimated_abun.tif');

%  compare endmember estimations
figure;
subplot_tight(2, 2, 1,[.08 .08]); plot(A); xlim([0 L]); title('Ground truth','fontsize',8); 
subplot_tight(2, 2, 2,[.08 .08]); plot(A_VCA); xlim([0 L]); title('A-VCA','fontsize',8);      
subplot_tight(2, 2, 3,[.08 .08]); plot(A_init); xlim([0 L]); title('A-init','fontsize',8);  
subplot_tight(2, 2, 4,[.08 .08]); plot(W); xlim([0 L]); title('A-SGSNMF','fontsize',8);  
drawnow; 
saveas(gcf,'../results/estimated_endm.tif');