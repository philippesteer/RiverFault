function Earthquake = fEarthquake(Mainshock,Aftershock)
%==========================================================================
global parFAULT parEQ;
%==========================================================================
% Concatenate
Earthquake.Mw = horzcat(Mainshock.Mw,vertcat(Aftershock.Mw)');
Earthquake.x  = horzcat(Mainshock.x,vertcat(Aftershock.x)');
Earthquake.y  = horzcat(Mainshock.y,vertcat(Aftershock.y)');
Earthquake.z  = horzcat(Mainshock.z,vertcat(Aftershock.z)');
Earthquake.g  = horzcat(Mainshock.g,vertcat(Aftershock.g)');
Earthquake.t  = horzcat(Mainshock.t,vertcat(Aftershock.t)');
Earthquake.dt = horzcat(Mainshock.dt,vertcat(Aftershock.dt)');

% Sort by time
[~,inew] = sort(Earthquake.t,'ascend');
Earthquake.Mw = Earthquake.Mw(inew);
Earthquake.x  = Earthquake.x(inew);
Earthquake.y  = Earthquake.y(inew);
Earthquake.z  = Earthquake.z(inew);
Earthquake.g  = Earthquake.g(inew);
Earthquake.t  = Earthquake.t(inew);
Earthquake.dt = Earthquake.dt(inew);

% % Periodic fault
% inew = find(Earthquake.x<0);            Earthquake.x(inew)=parFAULT.W+Earthquake.x(inew);
% inew = find(Earthquake.x>parFAULT.W);   Earthquake.x(inew)=Earthquake.x(inew)-parFAULT.W;
% inew = find(Earthquake.y<0);            Earthquake.y(inew)=parFAULT.L+Earthquake.y(inew);
% inew = find(Earthquake.y>parFAULT.L);   Earthquake.y(inew)=Earthquake.y(inew)-parFAULT.L;

% Non-periodic fault - Remove Earthquakes outside the fault bounding box
inew = find(Earthquake.x>0 & Earthquake.x<parFAULT.W & Earthquake.y>0 & Earthquake.y<parFAULT.L & Earthquake.t<=parEQ.T);
Earthquake.Mw = Earthquake.Mw(inew);
Earthquake.x  = Earthquake.x(inew);
Earthquake.y  = Earthquake.y(inew);
Earthquake.z  = Earthquake.z(inew);
Earthquake.g  = Earthquake.g(inew);
Earthquake.t  = Earthquake.t(inew);
Earthquake.dt = Earthquake.dt(inew);

% Save the initial coordinates of the earthquakes
Earthquake.xini = Earthquake.x;
Earthquake.yini = Earthquake.y;
Earthquake.zini = Earthquake.z;

disp('Earthquake frequency-magnitude computed ...')









