clear all
close all
clc

ts = 1e-4;
omega_nom = 2*pi*50;

% ddot x = 0;
% Transition matrix in continuous time domain
Aso = [0 1; 0 0];
% Transition matrix in discrete time domain (ts)
Adso = eye(2) + Aso * ts;
Cso = [1 0];
% poles to be placed in c.t.d.
p2place = [-1 -4] * omega_nom;
% poles to be placed in d.t.d.
p2placed = exp(p2place * ts);
% state observer feedback gains
Ldso = (acker(Adso',Cso',p2placed))';
komega = Ldso(2) / ts;
ktheta = Ldso(1) / ts;

%% example 2
p2place = [-5 -20] * omega_nom;
p2placed = exp(p2place * ts);
Ldso2 = (acker(Adso',Cso',p2placed))';
komega2 = Ldso2(2) / ts;
ktheta2 = Ldso2(1) / ts;