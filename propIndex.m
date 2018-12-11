function [propInd,indNd] = propIndex(Nd_tmp)


indNd = find(Nd_tmp > 0);

propInd = [];

for j = 1:length(indNd)
    if(Nd_tmp(indNd(j)) > 1)
        
        propInd = [propInd indNd(j)*ones(1,Nd_tmp(indNd(j)))];
        
    else
        
        propInd = [propInd indNd(j)];
        
    end
end

propInd = propInd';