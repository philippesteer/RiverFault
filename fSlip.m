function [Slip]=fSlip(Earthquake)
%==========================================================================
% Parameters
global parFAULT parRUP;

x=[0:parFAULT.res:parFAULT.W];y=[0:parFAULT.res:parFAULT.L];
[Slip.X,Slip.Y]=meshgrid(x,y);Slip.Z=zeros(size(Slip.X));Slip.D=zeros(size(Slip.X));Slip.N=zeros(size(Slip.X));

% if numel(Earthquake.x)<=numel(Slip.X)
%     for i=1:numel(Earthquake.x)
%         ind=find(Slip.X>=Earthquake.x(i)+Earthquake.W1(i) & Slip.X<=Earthquake.x(i)+Earthquake.W2(i) & Slip.Y>=Earthquake.y(i)+Earthquake.L1(i) & Slip.Y<=Earthquake.y(i)+Earthquake.L2(i));
%         Slip.D(ind)=Slip.D(ind)+Earthquake.D(i);
%     end
% else
%     for i=1:numel(Slip.X)
%         ind=find(Earthquake.x+Earthquake.W1<=Slip.X(i) & Earthquake.x+Earthquake.W2>=Slip.X(i) & Earthquake.y+Earthquake.L1<=Slip.Y(i) & Earthquake.y+Earthquake.L2>=Slip.Y(i));
%         Slip.D(i)=Slip.D(i)+sum(Earthquake.D(ind));
%     end
% end

% Periodic Fault
tempX10=Earthquake.x+Earthquake.W1;tempX11=Earthquake.x.*0-1e32;
tempX20=Earthquake.x+Earthquake.W2;tempX21=Earthquake.x.*0-1e32;
tempY10=Earthquake.y+Earthquake.L1;tempY11=Earthquake.y.*0-1e32;
tempY20=Earthquake.y+Earthquake.L2;tempY21=Earthquake.y.*0-1e32;
% Find out-bounded ruptures
ind=find(tempX10<0);            tempX21(ind)=parFAULT.W; tempX11(ind)=parFAULT.W+tempX10(ind);  tempX10(ind)=0;
ind=find(tempX20>parFAULT.W);   tempX11(ind)=0;          tempX21(ind)=tempX20(ind)-parFAULT.W;  tempX20(ind)=parFAULT.W;
ind=find(tempY10<0);            tempY21(ind)=parFAULT.L; tempY11(ind)=parFAULT.L+tempY10(ind);  tempY10(ind)=0;
ind=find(tempY20>parFAULT.L);   tempY11(ind)=0;          tempY21(ind)=tempY20(ind)-parFAULT.L;  tempY20(ind)=parFAULT.L;


% Determine cumulated slip on the 2D discretized fault plane
for i=1:numel(Slip.X)
    ind=find(  ( (tempX10<=Slip.X(i) & tempX20>=Slip.X(i)) | (tempX11<=Slip.X(i) & tempX21>=Slip.X(i)) )  &  ( (tempY10<=Slip.Y(i) & tempY20>=Slip.Y(i)) | (tempY11<=Slip.Y(i) & tempY21>=Slip.Y(i)) )   );
    Slip.D(i)=Slip.D(i)+sum(Earthquake.D(ind));
    Slip.N(i)=Slip.N(i)+numel(ind);
end
% Keep only Earthquakes that ruptured the surface
isur=find(Slip.X==0);
Slip.indsur=[];
for i=1:numel(isur)
    ind=find(  ( (tempX10<=Slip.X(isur(i)) & tempX20>=Slip.X(isur(i))) | (tempX11<=Slip.X(isur(i)) & tempX21>=Slip.X(isur(i))) )  &  ( (tempY10<=Slip.Y(isur(i)) & tempY20>=Slip.Y(isur(i))) | (tempY11<=Slip.Y(isur(i)) & tempY21>=Slip.Y(isur(i))) )   );
    Slip.indsur=[Slip.indsur ind];
end
Slip.indsur=unique(Slip.indsur);

% Keep only Earthquakes that ruptured the river
iriv=find(Slip.X==0 & Slip.Y==parFAULT.L/2);
Slip.indriv=[];
for i=1:numel(iriv)
    ind=find(  ( (tempX10<=Slip.X(iriv(i)) & tempX20>=Slip.X(iriv(i))) | (tempX11<=Slip.X(iriv(i)) & tempX21>=Slip.X(iriv(i))) )  &  ( (tempY10<=Slip.Y(iriv(i)) & tempY20>=Slip.Y(iriv(i))) | (tempY11<=Slip.Y(iriv(i)) & tempY21>=Slip.Y(iriv(i))) )   );
    Slip.indriv=[Slip.indriv ind];
end
Slip.indriv=unique(Slip.indriv);

% Keep only Earthquakes that ruptured the river - several rivers
dy_rivs=parFAULT.res.*2;
y_rivs=[0:dy_rivs:parFAULT.L];
for j=1:numel(y_rivs)
    irivs=find(Slip.X==0 & Slip.Y==y_rivs(j));
    Slip.indrivs{j}=[];
    for i=1:numel(irivs)
        ind=find(  ( (tempX10<=Slip.X(irivs(i)) & tempX20>=Slip.X(irivs(i))) | (tempX11<=Slip.X(irivs(i)) & tempX21>=Slip.X(irivs(i))) )  &  ( (tempY10<=Slip.Y(irivs(i)) & tempY20>=Slip.Y(irivs(i))) | (tempY11<=Slip.Y(irivs(i)) & tempY21>=Slip.Y(irivs(i))) )   );
        Slip.indrivs{j}=[Slip.indrivs{j} ind];
    end
    Slip.indrivs{j}=unique(Slip.indrivs{j});
end

