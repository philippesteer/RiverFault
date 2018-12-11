function [var] = fBASSlocal(MW,X,Y,Z,G,t,dt)
%==========================================================================
% Parameters
global parEQ parRUP;

%==========================================================================
% First generation of daughter EQ.
Nd{1}    = round(10.^(parEQ.b*(MW - parEQ.Dm - parEQ.mw_min)));
% Random vectors generation
rnd.Pm   = rand(Nd{1},1);
rnd.Pt   = rand(Nd{1},1);
rnd.Pr   = rand(Nd{1},1);
% Angular sampling - circular 
rnd.Pphi = 2.*pi.*rand(Nd{1},1);
% Magnitude, Occurrence Time, Coord and Generation number of first generation
mw_d{1}  = (log10(rnd.Pm) - parEQ.b * parEQ.mw_min)/-parEQ.b;
% Bounded magnitude -----------------------
mw_d{1}(mw_d{1}>parRUP.Mwmax)=parRUP.Mwmax;
%------------------------------------------
r_d{1}   = ((1 ./ rnd.Pr).^(1/(parEQ.q-1)) - 1).*parEQ.d.*10.^(0.5*MW);
x_d{1}   = X + r_d{1} .* cos(rnd.Pphi);
y_d{1}   = Y + r_d{1} .* sin(rnd.Pphi);
z_d{1}   = Z + zeros(Nd{1},1);
g_d{1}   = G + ones(Nd{1},1);
t_inc= parEQ.c.*((1 ./ rnd.Pt).^(1/(parEQ.p-1)) - 1);
t_d{1}   = t + t_inc;
dt_d{1}  = dt + t_inc;

% Variable initialization
var.Mw   = [mw_d{1}];
var.x    = [x_d{1}];
var.y    = [y_d{1}];
var.z    = [z_d{1}];
var.g    = [g_d{1}];
var.t    = [t_d{1}];
var.dt   = [dt_d{1}];

inc = 1;
while Nd{inc} > 0
    
    % Number of aftershocks of the inc-generation
    Nd_tmp{inc+1} = round(10.^(parEQ.b*(mw_d{inc} - parEQ.Dm - parEQ.mw_min)));
    Nd{inc+1}     = sum(Nd_tmp{inc+1});
    indNd         = propIndex(Nd_tmp{inc+1});
    
    if isempty(indNd)==0
        
        % Random vectors generation
        rnd.Pm      = rand(Nd{inc+1},1);
        rnd.Pt      = rand(Nd{inc+1},1);
        rnd.Pr      = rand(Nd{inc+1},1);
        % Angular sampling - circular 
        rnd.Pphi    = 2.*pi.*rand(Nd{inc+1},1);
        % Magnitude, Occurrence Time and Coord first generation
        mw_d{inc+1} = (log10(rnd.Pm) - parEQ.b * parEQ.mw_min)/-parEQ.b;
        % Bounded magnitude -----------------------
        mw_d{inc+1}(mw_d{inc+1}>parRUP.Mwmax)=parRUP.Mwmax;
        %------------------------------------------
        r_d{inc+1}  = ((1 ./ rnd.Pr).^(1/(parEQ.q-1)) - 1).* parEQ.d.*10.^(0.5*mw_d{inc+1});
        x_d{inc+1}  = x_d{inc}(indNd) + r_d{inc+1} .* cos(rnd.Pphi);
        y_d{inc+1}  = y_d{inc}(indNd) + r_d{inc+1} .* sin(rnd.Pphi);
        z_d{inc+1}  = z_d{inc}(indNd) + zeros(Nd{inc+1},1);
        g_d{inc+1}  = g_d{inc}(indNd) + ones(Nd{inc+1},1);
        t_inc=parEQ.c.*((1 ./ rnd.Pt).^(1/(parEQ.p-1)) - 1);
        t_d{inc+1}  = t_d{inc}(indNd) + t_inc;
        dt_d{inc+1} = dt_d{inc}(indNd)+ t_inc;

        % Variable update [push-back]
        var.Mw   = [var.Mw; mw_d{inc+1}];
        var.x    = [var.x; x_d{inc+1}];
        var.y    = [var.y; y_d{inc+1}];
        var.z    = [var.z; z_d{inc+1}];
        var.g    = [var.g; g_d{inc+1}];
        var.t    = [var.t; t_d{inc+1}];
        var.dt   = [var.dt;dt_d{inc+1}];
        
    end
       
    inc = inc + 1;
end



