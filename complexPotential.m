%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ASEN 5051
%                           Mini Project 1
%
% This script is meant to be used to calculate the complex potential
% associated with lifting flow over an airfoil (to be used in a Joukowski
% transfrom to model potential flow over an airfoil).
%
% Authors: Lucas Calvert, Duncan McGough
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
clear all;
close all;
clc;

%% Setup

syms   u_inf alpha R zeta zetaPrime Gamma b

%Enter stream function
Psi = u_inf*exp(-1i*alpha)*(zeta-zetaPrime) + u_inf*exp(1i*alpha)*(R^2)/(zeta-zetaPrime) + (Gamma/(2*pi*1i))*log((zeta-zetaPrime)/R);

%Take derivative
dPsi = diff(Psi,zeta)

J = zeta + b^2/zeta;

dJ = diff(J,zeta)


vel = dPsi/dJ

%Set derivative equal to zero and  solve for Gamma
%Gamma = solve(dPsi==0,Gamma)


