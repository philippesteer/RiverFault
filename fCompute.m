function [Mainshock,Earthquake,Slip]=fCompute()    

global  parFAULT parRUP parEQ ;

Mainshock  = fMainshock();
tic
Aftershock(1,numel(Mainshock)).Mw=[];
Aftershock(1,numel(Mainshock)).x=[];
Aftershock(1,numel(Mainshock)).y=[];
Aftershock(1,numel(Mainshock)).z=[];
Aftershock(1,numel(Mainshock)).g=[];
Aftershock(1,numel(Mainshock)).t=[];
Aftershock(1,numel(Mainshock)).dt=[];

for i = 1:length(Mainshock.Mw)
    Aftershock(i) = fBASS_local(Mainshock.Mw(i),Mainshock.x(i),Mainshock.y(i),Mainshock.z(i),Mainshock.g(i),Mainshock.t(i),Mainshock.dt(i));
end
toc
Earthquake = fEarthquake(Mainshock,Aftershock);

%% Earthquake Rupture & fault slip
[Earthquake] = fRupture(Earthquake);
[Slip]=fSlip(Earthquake);
