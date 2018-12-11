function [var] = frelocateonfault(var)
%==========================================================================
% Parameters
global parFAULT;

% rotate the Earthquake coordinates to respect the dip
x=var.x;y=var.y;z=var.z; % Stock coordinate
var.x = (x.*cosd(parFAULT.dip))-((z-parFAULT.ztop).*sind(parFAULT.dip)); 
var.y = y; 
var.z = (x.*sind(parFAULT.dip))+((z-parFAULT.ztop).*cosd(parFAULT.dip))+parFAULT.ztop; 
% rotate  Earthquake coordinates to respect the strike
x=var.x;y=var.y;z=var.z; % Stock coordinate
var.x = (x.*cosd(-parFAULT.strike))-(y.*sind(-parFAULT.strike)); 
var.y = (x.*sind(-parFAULT.strike))+(y.*cosd(-parFAULT.strike)); 
var.z = var.z;
% Translate coordinates
x=var.x;y=var.y;z=var.z; % Stock coordinate
var.x=x+parFAULT.xll;
var.y=y+parFAULT.yll;
var.z=z;