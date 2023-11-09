function [F, s] = flux_function(u, v_n, mi)
%eulerflux calculates the flux components via the governing equations from 
%the state u. Additionally outputs the ion sound speed c.
%   u = [n_n; n_e; n_e*v_e; P_e];
%   F = [n_n*v_n; n_e*v_e; n_e*v_e^2 + P_e/m_i; 5/2*P_e*v_e];
%   c = sqrt( (5/3*P_e) / (n_e*m_i) )

pe = u(4)/1.5;
v_e = u(3)./u(2);

c = sqrt(pe/ (u(2)*mi) ); %ion sound speed  
m_molar = 131.2930/1000; %kg/mol
c_n = sqrt( ((5/3)*8.3145*300)/m_molar); %J / kg
F = [u(1)*v_n,  u(3), ((u(3)*v_e) + (pe/mi)), (5/2)*pe*v_e]; %flux

s = [v_n+c_n, c+v_e, c+v_e, c+v_e]; %wave speeds

v_e = [v_n, v_e, v_e, v_e];
c = [c_n, c, c, c];

end