function fInput()

global  parFAULT parRUP parEQ ;


% Parameters defining fault extent and seismicity distribution
% Model parameters - Fault 
parFAULT.L      = 200000;
parFAULT.W      = 30000;
parFAULT.res    = 500;
% Distribution of seismicity along fault width
parFAULT.dist   = 'Normal'; % or 'Uniform'
parFAULT.mu     = parFAULT.W/2; % peak of seismicity
parFAULT.sigma  = parFAULT.W/5; % distribution half-width
parFAULT.chi=0.5; % fraction of seismic slip over total fault slip
% Compute the ratio between the Gaussian distribution and the normal
% distribution - parFAULT.ratio is the peak to peak ratio
if (strcmp(parFAULT.dist,'Normal'))
    pd = makedist('Normal','mu',parFAULT.mu,'sigma',parFAULT.sigma);pdN = truncate(pd,0,parFAULT.W);
    parFAULT.ratio=pdf(pdN, parFAULT.mu)./unifpdf(parFAULT.mu,0,parFAULT.W);
else
    parFAULT.ratio=1;
end
%% 
% Parameters defining ruptures and the used scaling laws

% Model parameters - Earthquake rupture [Leonard, 2010]
% Reverse fault
parRUP.mu       = 30e9;                            % Shear modulus (Pa)
parRUP.Beta     = 2/3;                             % parameter Leonard or Beta=1
parRUP.C1       = 17.5;                            % parameter Leonard
parRUP.C2       = 3.8e-5;                          % parameter Leonard
% Define maximum and minimum earthquake magnitude
parRUP.Lmax     = parFAULT.L;
parRUP.Wmax     = parFAULT.W;
parRUP.MoLmax   = parRUP.mu*parRUP.C1^(3/2)*parRUP.C2*parFAULT.L.^(3/2*(1+parRUP.Beta));
parRUP.MoWmax   = parRUP.mu*parRUP.C1^(-3/(2*parRUP.Beta))*parRUP.C2*parFAULT.W.^(3/2*(1+1/parRUP.Beta));
parRUP.Momax    = min(parRUP.MoLmax,parRUP.MoWmax);
parRUP.Mwmax    = 2/3*log10(parRUP.Momax)-6.07;
parRUP.Lmin     = parFAULT.res;
parRUP.Wmin     = parFAULT.res;
parRUP.MoLmin   = parRUP.mu*parRUP.C1^(3/2)*parRUP.C2*parRUP.Lmin.^(3/2*(1+parRUP.Beta));
parRUP.MoWmin   = parRUP.mu*parRUP.C1^(-3/(2*parRUP.Beta))*parRUP.C2*parRUP.Wmin.^(3/2*(1+1/parRUP.Beta));
parRUP.Momin    = max(parRUP.MoLmin,parRUP.MoWmin);
parRUP.Mwmin    = 2/3*log10(parRUP.Momin)-6.07;
%% 
% Parametersused in the earthquake generator model

% Model parameters - EarthQuake generator [Turcotte et al., 2007]
parEQ.type      = 'Gutenberg-Richter';
parEQ.T         = 2e3.*365;                                                % Total time duration (days)
parEQ.mw_min    = parRUP.Mwmin;                                            % Minimum Magnitude (deterministic)
parEQ.mw_max    = parRUP.Mwmax;                                            % Maximum Magnitude (deterministic)
parEQ.mw_dm     = 0.02;                                                    % Magnitude steps
parEQ.mw        = parEQ.mw_min:parEQ.mw_dm:parEQ.mw_max;
parEQ.b         = 1;                                                       % GR param for Mainshocks and Aftershocks
parEQ.Dm        = 1.25;                                                    % GR-Bath param
parEQ.c         = 0.1;                                                     % temporal Omori param
parEQ.p         = 1.25;                                                    % temporal Omori param
parEQ.d         = 4.0;                                                     % spatial Omori param
parEQ.q         = 1.35;                                                    % spatial Omori param
if (strcmp(parEQ.type,'Periodic'))
    parEQ.period    = 330*365;                                                 % Periodicity of the Main shocks
    parEQ.Mwevent   = 8;                                                       % Magnitude of the reccurent event
    parEQ.Nmain     = floor(parEQ.T./parEQ.period)+1;                          % Number of mainshocks
elseif (strcmp(parEQ.type,'Gutenberg-Richter'))
    parEQ.srate     = 0.1./parFAULT.ratio;                                     % Number of mainshock per day
    parEQ.Nmain     = parEQ.T*parEQ.srate;                                     % Number of mainshocks
    parEQ.a         = log10(parEQ.T.*parEQ.srate) + parEQ.b*parEQ.mw_min;      % GR param for Mainshocks - for Mw>=Mwmin
end

disp('Input parameters - done')
