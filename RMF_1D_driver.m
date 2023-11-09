close all; clear all; clc;

%Thruster parameters
m_dot = 272/11.12e6; %kg/s
l = 0.11; %m
r_0 = 0.07; %m %throat radius
theta = 46; %deg cone half angle (35)
Br = 43e-4; %Tesla
f_rmf = 413.2; %kHz

[eta,Thrust,Isp,ne,Te,nn] = RMF_1D_OG(m_dot,f_rmf,Br,l,r_0,theta);
