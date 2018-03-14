% Show segmentation results
% Usage: Img_seg = showsegresults(Y,labels);
%==========================================================================
function Img_seg = showsegresults(Y,labels)
[p,q,L] = size(Y);
Img=zeros(p,q,3); lim = [round(L/4),round(L/2),round(3*L/4)];
for itr= 1:3
    Sp = Y(:,:,lim(itr));
    X_r = sort(reshape(Sp,p*q,1)); 
    a=X_r(round(0.02*p*q));
    b=X_r(round(0.98*p*q));
    for i=1:p
        for j=1:q
            if Sp(i,j) < a
                Img(i,j,itr) = a;
            else
                if Sp(i,j)< b
                Img(i,j,itr) = Sp(i,j);
                else
                Img(i,j,itr) = b;
                end
            end
        end    
    end
    Img(:,:,itr)  =(Sp -a).*(255/(b-a));
end
Img=uint8(Img);
Img_seg = drawregionboundaries(labels, Img, [255 0 0]);
end



