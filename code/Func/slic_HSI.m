% Preprocessing of SGSNMF:
% Simple Linear Iterative Clustering SuperPixels (SLIC) for HSIs
%
% Usage:   [l, Am, S ,C, Cj] = slic_HSI(im, k, m, seRadius, nItr)
%==========================================================================
% Input:      im - 3D hyperspectral image to be segmented.
%              k - Number of desired superpixels.
%              m - Weighting factor between spectral and spatial differences. {0.5}
%       seRadius - Regions morphologically smaller than this are merged with adjacent regions. {1}
%           nItr - the number of iterations
%
% Output:      l - Labeled image of superpixels. Labels range from 1 to k.
%             Am - Adjacency matrix of segments.  Am(i, j) indicates whether
%                  segments labeled i and j are connected/adjacent
%              S - the averange size of superpixels
%              C - cluster center of each superpixel
%             Cj - the confidence index
%==========================================================================
% 
% Copyright (C) 2018, Xinyu Wang (wangxinyu@whu.edu.cn)
%                     Yanfei Zhong (zhongyanfei@whu.edu.cn)
%                     Wuhan university
%                     All rights reserved.
%
% Reference:
%
% X. Wang, Y. Zhong, L. Zhang, and Y. Xu, ¡°Spatial Group Sparsity
% Regularized Nonnegative Matrix Factorization for Hyperspectral Unmixing,¡±
% IEEE Transactions on Geoscience and Remote Sensing, vol. 55, no. 11, pp.
% 6287-6304, NOV 2017, 2017.
%
% R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and S. Susstrunk. "SLIC
% Superpixels Compared to State-of-the-Art Superpixel Methods"  PAMI. Vol
% 34 No 11. November 2012. pp 2274-2281. 
%
% Last Modified:
% 27 Feb, 2018

function [l, Am, S ,C ,Cj] = slic_HSI(im, k, m, seRadius, nItr)

%% Initialization
    if ~exist('seRadius','var')   || isempty(seRadius),     seRadius = 1;     end
    if ~exist('nItr','var')   || isempty(nItr),     nItr = 10;     end
    [rows, cols, L] = size(im);

    S = sqrt(rows*cols / (k * sqrt(3)/2));
    countCols = round(cols/S - 0.5);
    S = cols/(countCols + 0.5); 
    countRows = round(rows/(sqrt(3)/2*S));
    vSpacing = rows/countRows;

    % Recompute the number of superpixels k
    k = countRows * countCols;
    
    % Allocate memory and initialise clusters, labels and distances.
    C = zeros(L+3,k);        % Cluster centre data  1:L is mean spectral value,
                             % L+1: L+2 is row, col of centre, 
                             % L+3 is the number of pixels in each clustering
    l = -ones(rows, cols);   % Pixel labels.
    d = inf(rows, cols);     % Pixel distances from cluster centres.
    
    % Initialise clusters on a hexagonal grid
    kk = 1;
    r = vSpacing/2;
    for ri = 1:countRows
        if mod(ri,2), c = S/2; else, c = S;  end
        for ci = 1:countCols
            cc = round(c); rr = round(r);
            C(1:L+2, kk) = [squeeze(im(rr,cc,:)); cc; rr];
            c = c+S;
            kk = kk+1;
        end
        r = r+vSpacing;
    end
    S = round(S);
    
 %  Upadate cluster centers   
    for n = 1:nItr
       for kk = 1:k
           % Get subimage around cluster
           rmin = max(C(L+2,kk)-S, 1);   rmax = min(C(L+2,kk)+S, rows); 
           cmin = max(C(L+1,kk)-S, 1);   cmax = min(C(L+1,kk)+S, cols); 
           subim = im(rmin:rmax, cmin:cmax, :);  
           assert(numel(subim) > 0);
           % Calculate distance
           D = dist(C(:, kk), subim, rmin, cmin, S, m);
           % If any pixel distance from the cluster centre is less than its
           % previous value update its distance and label
           subd =  d(rmin:rmax, cmin:cmax);
           subl =  l(rmin:rmax, cmin:cmax);
           updateMask = D < subd;
           subd(updateMask) = D(updateMask);
           subl(updateMask) = kk;
           d(rmin:rmax, cmin:cmax) = subd;
           l(rmin:rmax, cmin:cmax) = subl;           
       end
       % Update cluster centres with mean spectral values
       C(:) = 0;
       for r = 1:rows
           for c = 1:cols
              spectra = reshape(im(r,c,:),L,1);
              tmp = [spectra; c; r; 1];
              C(:, l(r,c)) = C(:, l(r,c)) + tmp;
           end
       end
       % Divide by number of pixels in each superpixel to get mean spectrum
       for kk = 1:k 
           C(1:L,kk) = C(1:L,kk)/C(L+3,kk); 
           C(L+1:L+2,kk) = round(C(L+1:L+2,kk)/C(L+3,kk)); 
       end
    end
%% Cleanup small orphaned regions.     
    % The cleaned up regions are assigned to the nearest cluster.
    if seRadius
        [l, Am] = mcleanupregions(l, seRadius);
    else
        l = makeregionsdistinct(l);
        l = renumberregions(l);
        Am = regionadjacency(l);    
    end
%% recalculate the center
    N = length(Am);
    C = zeros(L+3,N); 
    for r = 1:rows
       for c = 1:cols
          spec = reshape(im(r,c,:),L,1);
          tmp = [spec; c; r; 1];
          C(:, l(r,c)) = C(:, l(r,c)) + tmp;
       end
    end
    % Divide by number of pixels in each superpixel to get mean values
    for kk = 1:N
       C(1:L,kk) = C(1:L,kk)/C(L+3,kk); 
       C(L+1:L+2,kk) = round(C(L+1:L+2,kk)/C(L+3,kk)); 
    end
    % Calculate the confidence index
    d_c = zeros(rows,cols); % Spectral distance
    d_s = zeros(rows,cols); % Spatial distance
    for r = 1:rows
        for c = 1:cols
            lab = l(r,c); 
            Cspec = C(1:L,lab);
            spec= reshape(im(r,c,:),L,1);
            d_c(r,c) = acos(min(max(dot(spec,Cspec)/(norm(spec)*norm(Cspec)),0.00000),1.000000));
            d_s(r,c) = (r - C(L+2,lab)).^2+(c - C(L+1,lab)).^2;
        end
    end
    Cj = 1./(sqrt(d_c + d_s/S.^2*m)+eps);
end
       
    
%%%%%%%%%%%%%%%%%%%%%%  dist  %%%%%%%%%%%%%%%%%%%%%%
%
% Usage:  D = dist(C, im, r1, c1, S, m)
% 
% Arguments:   C - Cluster center
%             im - sub-image surrounding cluster center
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing
%              m - weighting factor between spectral and spatial distances.
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster center
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = arccos(ims*imc/sqrt(sum(ims.^2,2)*sum(imc.^2)))  % Spectral distance
% ds = sqrt(dx^2 + dy^2)                                % Spatial distance

function D = dist(C, im, r1, c1, S, m)
    % Squared spatial distance
    [rows, cols, L] = size(im);
    [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
    x = x-C(L+1);  % x and y dist from cluster centre
    y = y-C(L+2);
    ds2 = x.^2 + y.^2;
    % Squared colour difference
    ims = reshape(im,rows*cols,L);
    imc = C(1:L);
    dcos = (ims*imc)./sqrt(sum(ims.^2,2)*sum(imc.^2));
    dc2 = acos((min(max(reshape(dcos,rows,cols),0.0000),1.0000)));
    D = sqrt(dc2 + ds2/S^2*m);
end
