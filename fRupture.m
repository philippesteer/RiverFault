function [Earthquake] = fRupture(Earthquake)
%==========================================================================
% Parameters
global parFAULT parRUP;

% Determine the geometrical properties of earthquakes (Leonard 2010)
Earthquake.Mo = 10.^(3/2.*Earthquake.Mw+9.105); 
    
% Leonard (2010)- Reverse fault - General case
Earthquake.L  = (Earthquake.Mo./(parRUP.mu*parRUP.C1^(3/2)*parRUP.C2)).^(2/(3*(1+parRUP.Beta)));
Earthquake.W  = parRUP.C1.*Earthquake.L.^parRUP.Beta;
Earthquake.D  = parRUP.C2*parRUP.C1^(1/2)*Earthquake.L.^(0.5*(1+parRUP.Beta));

% Randomly distribute rupture around the hypocenter
Wdown=rand(1,numel(Earthquake.W)).*Earthquake.W; Wup=Earthquake.W-Wdown;
Lleft=rand(1,numel(Earthquake.L)).*Earthquake.L; Lright=Earthquake.L-Lleft;

% % Readjust rupture patch location if it exceeds fault dimension
% ind=find(Earthquake.x-Wdown<0);Wdown(ind)=Earthquake.x(ind);Wup=Earthquake.W-Wdown;
% ind=find(Earthquake.x+Wup>parFAULT.W);Wup(ind)=parFAULT.W-Earthquake.x(ind);Wdown=Earthquake.W-Wup;
% ind=find(Earthquake.y-Lleft<0);Lleft(ind)=Earthquake.y(ind);Lright=Earthquake.L-Lleft;
% ind=find(Earthquake.y+Lright>parFAULT.L);Lright(ind)=parFAULT.L-Earthquake.y(ind);Lleft=Earthquake.L-Lright;

% Stock information
Earthquake.W1=-Wdown;Earthquake.W2=Wup;
Earthquake.L1=-Lleft;Earthquake.L2=Lright;

