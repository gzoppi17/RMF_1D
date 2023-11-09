function [F_F] = russonov(u_l,u_j,u_r, v_n, mi)
%computes net russonov flux for current cell state u_j, and neighbor
%states u_l, and u_r.

%compute fluxes from left and right cells
[F_l,s_l] = flux_function(u_l, v_n, mi);
[F, s] = flux_function(u_j, v_n, mi);
[F_r,s_r] = flux_function(u_r, v_n, mi);

s_max_l = max([s_l; s]);
s_max_r = max([s; s_r]);

a = max(max([s_max_l; s_max_r]));

F_in = 0.5*(F_l+F)-0.5*s_max_l.*(u_j-u_l);
F_out = 0.5*(F+F_r)-0.5*s_max_r.*(u_r-u_j);

F_F = F_in - F_out;

end