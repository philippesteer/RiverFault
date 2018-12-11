function Mainshock = fMainshock()
%==========================================================================
% Generate mainshock sequence
%==========================================================================
% Load relevant parameters
global parEQ parFAULT;

if (strcmp(parEQ.type,'Periodic'))
 
    % Mainshock magnitude
    Mainshock.Mw=ones(1,parEQ.Nmain).*parEQ.Mwevent;
    
    % Mainshock time
    Mainshock.t   = 0:parEQ.period:parEQ.T;
          
elseif (strcmp(parEQ.type,'Gutenberg-Richter'))
    
    % Mainshock magnitude
    Mainshock.Mw  = [];
    for i = 1:length(parEQ.mw)
        N            = round(10^(parEQ.a-parEQ.b*parEQ.mw(i)))-round(10^(parEQ.a-parEQ.b*(parEQ.mw(i)+parEQ.mw_dm)));
        Mainshock.Mw = [Mainshock.Mw parEQ.mw(i) + parEQ.mw_dm.*(rand(1,N)-0.5)];
    end
    
    % Shuffle vector
    Mainshock.Mw  = Mainshock.Mw(randperm(length(Mainshock.Mw)));
    
    % Mainshock time
    Mainshock.t   = parEQ.T.*rand(1,length(Mainshock.Mw));

end    
% Mainshock Coordinates epicenter

dx            = 0;
if (strcmp(parFAULT.dist,'Normal'))
    pd = makedist('Normal','mu',parFAULT.mu,'sigma',parFAULT.sigma);pd = truncate(pd,0,parFAULT.W);
    Mainshock.x   = random(pd,[1,length(Mainshock.Mw)]);
elseif (strcmp(parFAULT.dist,'Uniform'))    
    pd = makedist('Uniform','lower',0,'upper',parFAULT.W);%pd = truncate(pd,0,parFAULT.W);
    Mainshock.x   = random(pd,[1,length(Mainshock.Mw)]);   
    %Mainshock.x   = round(dx + ((parFAULT.W-dx) - dx).*rand(1,length(Mainshock.Mw)));
end
pd = makedist('Uniform','lower',0,'upper',parFAULT.L);
Mainshock.y   = random(pd,[1,length(Mainshock.Mw)]);   
Mainshock.z   = zeros(1,length(Mainshock.Mw));  

% Initialize - Mainshock aftershock generation number (by definition 0)
Mainshock.g   = zeros(1,length(Mainshock.Mw));

% Initialize -  Mainshock dt
Mainshock.dt  = zeros(1,length(Mainshock.Mw));
