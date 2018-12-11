% ## RMSE
% Usage: rmse = rmse(S_GT,S_E,Ab_ord,S_E_ord)
%======================================================================

function rmse = rmse(S_GT,S_E,GT_ord,E_ord)

[P,N]= size(S_E);
rmse= zeros(3,P+1);
rmse(1,:) = GT_ord;
rmse(2,:) = E_ord;
for i=1:P
    rmse(3,i) = sqrt(sum((S_GT(i,:)- S_E(rmse(2,i),:)).^2)/N);
end
rmse(3,P+1)= sqrt(sum(rmse(3,1:P).^2));

fprintf([repmat('%8i',1,P+1),'\n'],rmse(1:2,:)');
fprintf([repmat('%8.4f',1,P+1),'\n'],rmse(3,:)');
end