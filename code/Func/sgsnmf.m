% Spatial Group Sparsity Regularized Nonnegative Matrix Factorization for
% Hyperspectral Unmixing
%
% Usage: [W,H] = sgsnmf(para,seg);
%==========================================================================
% Input:  para.X (V)  -  L*N hyperspectrtal data;
%         para.W (W)  -  L*M initial endmember matrix;
%         para.H (H)  -  M*N initial abundance matrix; 
%        para.lambda  -  regularization parameter; default{0.3} 
%           para.tol  -  tolerance for a relative stopping condition; default{0.05}
%       para.maxiter  -  Stop condition 1: limit of iterations; default{100}
%     para.timelimit  -  Stop condition 2: limit of time;
%       para.verbose  -  1- display current information in SGSNMF, and 0 - otherwise  
%    para.print_iter  -  iterations between print on screen
%
%              seg.P  -  the number of superpixels
%            seg.X_c  -  the averange spectra of superpixels 
%             seg.Cj  -  the confidence index
%         seg.labels  -  the labels of pixels
%         
% Output: W  -  L*M estimated endmember matrix obtained by SGSNMF
%         H  -  M*N estimated abundance matrix obtained by SGSNMF
%==========================================================================
%
% Copyright (C) 2018, Xinyu Wang (wangxinyu@whu.edu.cn)
%                     Yanfei Zhong (zhongyanfei@whu.edu.cn)
%                     Wuhan university
%                     All rights reserved.
%
%
% References:
% X. Wang, Y. Zhong, L. Zhang, and Y. Xu, ¡°Spatial Group Sparsity
% Regularized Nonnegative Matrix Factorization for Hyperspectral Unmixing,¡±
% IEEE Transactions on Geoscience and Remote Sensing, vol. 55, no. 11, pp.
% 6287-6304, NOV 2017, 2017.

%
% Last Modified:
% 27 Feb, 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W,H] = sgsnmf(para,seg)
    
    % calculate initial gradient
    V = para.X; W = para.W; H = para.H; initt = cputime;
    gradW = W*(H*H') - V*H'; gradH = (W'*W)*H - W'*V;
    initgrad = norm([gradW; gradH'],'fro');
    fprintf('Init gradient norm %f\n', initgrad); 
    
    %tolerance for a relative stopping condition
    tolW = max(0.001,para.tol)*initgrad; 
    tolH = ones(seg.P,1).*tolW; 
    
    objhistory =  sum(sum((V-W*H).^2));
    objhistory = [objhistory 0];
    inc = 0;
    inc0 = 0;
    for iter=1:para.maxiter 
      % stopping condition
        % condition 1 time
        if cputime-initt > para.timelimit, break; end
        % condition 2 
        if objhistory(end-1)-objhistory(end)>0.0001
            inc = 0;
        else
            disp('inc');
            inc = inc+1;  
            inc0 = inc0+1;
        end
        if iter < 5, inc = 0; end 
        if inc >= 5 && inc0 >= 20, break; end     
        
       %  Update Wp
%         Wpmatrix =  hyperFcls(seg.X_c,W);
        Wpmatrix = fcls(W,seg.X_c);
        Wpmatrix = 1./(para.M.^2*Wpmatrix + 1);
        
       %  Update H    
        tW = ASC(W,15); tV = ASC(V,15); % ASC
        for P = 1 : seg.P
            labmask = (seg.labels == P); Cj = seg.Cj(labmask)';
            Wp = diag(Wpmatrix(:,P));
            [H(:,labmask),gradH(:,labmask),iterH] = subprobH(tV(:,labmask),tW,H(:,labmask),Wp,Cj,para.lambda,tolH(P),100);
            if iterH == 1, tolH(P) = 0.1 * tolH(P); end
        end
        
        %  Update W
        [W,gradW,iterW] = subprobW(V',H',W',tolW, 100); W = W'; gradW = gradW';
        if iterW==1
            tolW = 0.1 * tolW;
        end
        
       projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)],'fro');
       if para.verbose && rem(iter,para.print_iter)==0
            fprintf('\nIter = %d proj-grad norm %f\n', iter, projnorm);
       end
       objhistory = [objhistory sum(sum((V-W*H).^2))];
    end
    fprintf('\nIter = %d final proj-grad norm %f\n', iter, projnorm);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MA] = ASC(M,thresh)
    [L,P] = size(M);
    MA = thresh*ones(L+1,P);
    MA(1:L,:) = M;
end

function [H,grad,iter] = subprobW(V,W,Hinit,tol,maxiter)
    
    % H, grad: output solution and gradient
    % iter: #iterations used
    % V, W: constant matrices
    % Hinit: initial solution
    % tol: stopping tolerance
    % maxiter: limit of iterations
    
    H = Hinit; WtV = W'*V; WtW = W'*W; 
    alpha = 1; beta = 0.1;
    for iter=1:maxiter  
        grad = WtW*H - WtV;
        projgrad = norm(grad(grad < 0 | H >0));
        if projgrad < tol
            break;
        end
    % search step size 
        for inner_iter=1:20
            Hn = max(H - alpha*grad, 0); d = Hn-H;
            gradd=sum(sum(grad.*d)); dQd = sum(sum((WtW*d).*d));
            suff_decr = 0.99*gradd + 0.5*dQd < 0;
            if inner_iter==1, decr_alpha = ~suff_decr; Hp = H; end
            if decr_alpha  
                if suff_decr, H = Hn; break;
                else
                    alpha = alpha * beta;
                end
            else
                if ~suff_decr | Hp == Hn, H = Hp; break;
                else
                    alpha = alpha/beta; Hp = Hn;
                end
            end
        end
    end
    if iter==maxiter,  fprintf('Max iter in subprobW\n'); end
end

function [H,grad,iter] = subprobH(V,W,Hinit,Wp,Cj,lambda,tol,maxiter)

% H, grad: output solution and gradient
% iter: #iterations used
% V, W: constant matrices
% Hinit: initial solution
% tol: stopping tolerance
% maxiter: limit of iterations

    H = Hinit; WtV = W'*V; WtW = W'*W; WpW = Wp'*Wp;
    alpha = 1; beta = 0.1;
    for iter=1:maxiter
        par = Cj./sqrt(sum((Wp*H).^2,1));
        Bp = diag(par);
        grad = WtW*H - WtV + lambda* WpW *H *Bp;
        projgrad = norm(grad(grad < 0 | H >0));
        if projgrad < tol,  break;  end
        % search step size 
        for inner_iter=1:20 
            Hn = max(H - alpha*grad, 0); d = Hn-H;
            gradd=sum(sum(grad.*d)); dQd = sum(sum((WtW*d).*d));
            suff_decr = 0.99*gradd + 0.5*dQd < 0;
            if inner_iter==1, decr_alpha = ~suff_decr; Hp = H; end
            if decr_alpha  
                if suff_decr, H = Hn; break;
                else
                    alpha = alpha * beta;
                end
            else
                if ~suff_decr | Hp == Hn, H = Hp; break;
                else
                    alpha = alpha/beta; Hp = Hn;
                end
            end
       end
    end
%     if iter==maxiter, fprintf('Max iter in subprobH\n'); end
end