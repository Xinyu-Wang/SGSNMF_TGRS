% ## SAD
% Usage: [SAD, table_acos] = sam(A_GT, A_E)
%======================================================================

function [SAD, table_acos] = sam(A_GT, A_E)
[~,N_e] = size(A_E);
[~,N_t] = size(A_GT);

table_cos = (A_E'*A_GT)./sqrt(sum(A_E.^2)'*sum(A_GT.^2));
table_acos = acos(min(max(table_cos,0.0000),1.0000));
table3 = table_acos;
if N_e < N_t
    SAD = zeros(4,N_t+1);
	SAD(1,:)=[1:N_t 0];
    for j=1:N_e   
        sam_j = min(table3(:));
        [m,n]= find(table3 == sam_j);
        SAD(2,n(1)) = m(1);
        SAD(3,n(1)) = sam_j;
        SAD(4,n(1)) = sam_j*180/pi;
        table3(m(1),:)= 100;
        table3(:,n(1))= 100;
    end
    SAD(3,N_t+1)=sum(SAD(3,1:N_t))/N_e;
    SAD(4,N_t+1)=sum(SAD(4,1:N_t))/N_e;
	format_1 = [repmat('%8i',1,N_t+1),'\n'];
	format_2 = [repmat('%8.4f',1,N_t+1),'\n'];
	fprintf(format_1,SAD(1:2,:)');
	fprintf(format_2,SAD(3:4,:)');
else 
	SAD = zeros(5,N_t+1);
    for j=1:N_t   
        SAD(1,j)=j;
        sam_j = min(table3(:));
        [m,n]= find(table3 == sam_j);
        SAD(2,n(1)) = m(1);
        SAD(3,n(1)) = sam_j;
        SAD(4,n(1)) = sam_j*180/pi;
        table3(m(1),:)= 100;
        table3(:,n(1))= 100;
    end
    SAD(3,N_t+1)=sum(SAD(3,1:N_t))/N_t;
    SAD(4,N_t+1)=sum(SAD(4,1:N_t))/N_t;
    SAD(5,1:N_t)= sqrt(sqrt(sum(A_E(:,SAD(2,1:N_t)).^2)./sum(A_GT(:,SAD(1,1:N_t)).^2)));
    format_1 = [repmat('%8i',1,N_t+1),'\n'];
    format_2 = [repmat('%8.4f',1,N_t+1),'\n'];
    fprintf(format_1,SAD(1:2,:)');
    fprintf(format_2,SAD(3:5,:)');
end
end