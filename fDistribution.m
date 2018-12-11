function [Distrib]=fDistribution(Earthquake)
%==========================================================================
% Parameters
global parFAULT parRUP;

x=[parFAULT.res:parFAULT.res:parFAULT.W-parFAULT.res];y=[parFAULT.res:parFAULT.res:parFAULT.L-parFAULT.res];
[Distrib.X,Distrib.Y]=meshgrid(x,y);Distrib.Z=zeros(size(Distrib.X));Distrib.N=zeros(size(Distrib.X));
for i=1:numel(Distrib.X)
    ind=find(Earthquake.x>=Distrib.X(i)-parFAULT.res & Earthquake.x<Distrib.X(i)+parFAULT.res & Earthquake.y>=Distrib.Y(i)-parFAULT.res & Earthquake.y<Distrib.Y(i)+parFAULT.res);
    Distrib.N(i)=numel(ind);
end