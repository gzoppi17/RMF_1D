function [F_F, a, c, dpe] = HLLE_flux(u_l,u_j,u_r, v_n, mi, i)
%computes net russonov flux for current cell state u_j, and neighbor
%states u_l, and u_r.

%compute fluxes from left and right cells

[F_l,s_l,c_l, v_e_l] = flux_function(u_l,v_n, mi);
[F, s, c, v_e_j] = flux_function(u_j, v_n, mi);
[F_r,s_r, c_r, v_e_r] = flux_function(u_r, v_n, mi);


s_min_ll = min([[0,0,0,0]; v_e_l-c_l]);
s_min_lr = min([[0,0,0,0]; v_e_j-c]);

s_max_ll = max([[0,0,0,0]; v_e_l+c_l]);
s_max_lr = max([[0,0,0,0]; v_e_j+c]);

s_min_l = min([s_min_ll; s_min_lr]);
s_max_l = max([s_max_ll; s_max_lr]);

s_min_rl = min([[0,0,0,0]; v_e_j-c]);
s_min_rr = min([[0,0,0,0]; v_e_r-c_r]);

s_max_rl = max([[0,0,0,0]; v_e_j+c]);
s_max_rr = max([[0,0,0,0]; v_e_r+c_r]);

s_min_r = min([s_min_rl; s_min_rr]);
s_max_r = max([s_max_rl; s_max_rr]);

%a_l = abs(v_e_l)+c_l;
a = abs(v_e_j)+c;
a = max(a);
%a_r = abs(v_e_r)+c_r;
%a = max([a_l, a_r]);

F_in = 0.5.*(F_l+F)-0.5.*((s_max_l+s_min_l) ./ (s_max_l-s_min_l)).*(F-F_l)+ ((s_max_l.*s_min_l) ./ (s_max_l-s_min_l)).*(u_j-u_l);
F_out = 0.5.*(F+F_r)-0.5.*((s_max_r+s_min_r) ./ (s_max_r-s_min_r)).*(F_r-F)+ ((s_max_r.*s_min_r) ./ (s_max_r-s_min_r)).*(u_r-u_j);

F_F = F_in - F_out;

c = c(2);

dpe = u_r(4) - u_l(4);

end