function [ y, n ,Cn ] = addNoise( syntheticImage, noise_type ,SNR, eta, verbose)
%ADDNOISE Summary of this function goes here
%   Detailed explanation goes here
[m,n,L] = size(syntheticImage);
N = m*n;
x = reshape(syntheticImage,N,L)';
% adding noise
if verbose,fprintf(1,'Signal-to-noise ratio set to: %d dB\n', SNR);end
if isinf(SNR), n=zeros(L,N); Cn=zeros(L);      % no noise case
else
   if verbose, fprintf(1,'Noise type: %s\n',noise_type);end
   varianc = sum(x(:).^2)/10^(SNR/10) /L/N ;
   if strcmp(noise_type,'additive')              % color & white noise case
      % note that eta is 1/eta in equation (4) 
      quad_term = exp(-((1:L)-L/2).^2*eta^2/2);
      varianc_v = varianc*L*quad_term/sum(quad_term);
      Cn = diag(varianc_v); 
      n = sqrtm(Cn)*randn([L N]);
   end
   if strcmp(noise_type,'poisson')
      factor = varianc/mean(x(:)); 
      n = sqrt(factor) * sqrt(x) .* randn([L N]);
      Cn = diag(factor * mean(x,2)); %n*n'/N;
   end
end
y = x + n; % hyperspectral observations: equation(1)
return
end

