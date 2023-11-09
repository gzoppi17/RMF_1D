function [eta,Thrust,Isp,ne,Te,nn] = RMF_1D_OG(m_dot,f_rmf,Br,l,r_0,theta)
%set global valiables

global v_n mi f_c w

%Constants
e = 1.6022e-19; %c
me = 9.109383632e-31; %kg
mi = 2.18017e-25; %kg
mu0 = 1.2566e-6; %H/m
kb = 1.381e-23; %J/K
v_n = 1/4*sqrt((8*kb*300)/(pi*mi)); %m/s

Rppu = .200; %ohms
f_c = 0.01; %magnetic confinement (low is more confined

w = 2*pi*f_rmf*1000; %rad/s;

%Time stepping
dt_max = 1e-7;
dt = dt_max;
    
%% INITIALIZATION
Ncells = 25;

% %Set up grid

dx = l/(Ncells);
x = -dx:dx:l+dx; %cell boundary x-values
x_cell_edges = x(2:end-1);
x_cell = movmean(x_cell_edges,[0 1]);
x_cell = x_cell(1:Ncells);

%quasi-1D area
r_cell_edges = r_0 + x_cell_edges.*tand(theta);
r_cell = movmean(r_cell_edges,[0 1]);
r_cell = r_cell(1:end-1);

r = r_cell;
x = x_cell;

A = pi.*r.*r;
dAdx = 2*r*pi*tand(theta);

r_bar = 1/l*trapz(x,r); %TATE calculated value

%Set infow bc
n_n_0 = m_dot/(pi*r_0^2*v_n*mi);
n_e_0 = 1e16; 
v_e_0 = v_n;
P_e_0 = 2*e*n_e_0;

u_in = [n_n_0, n_e_0, n_e_0*v_e_0, 1.5*P_e_0]; 

%Set IC
u = u_in.*ones(Ncells,4);

u_next = zeros(Ncells,4); %cell update vector
R = zeros(Ncells,4); %residual vector
P = zeros(Ncells,1); %list of RMF power per unit length
F_F = zeros(Ncells,4); %flux vector
S_area = zeros(Ncells,4); %source vector
S = zeros(Ncells,4); %area expansion vector


%% MAIN ITERATION LOOP
count_max = 1000; %counter parameters for periodic state output
count = 0;
SumR = 10;
N = 0;

while sqrt(SumR) > 1e-6
    %% SINGLE TIME STEP 
    R = R.*0; %initialize residuals
    u_next = u_next.*0; %initralize update state
    count = count+1; %add value to counter
    N = N+1;

    for i = 1:Ncells %Loop over cells
        u_j = u(i,:); %state on current cell

        %Get neighbor states
        if i == 1 %if first cell 
            
            u_l = u_in; %use inflow state

        else
            u_l = u(i-1,:); %use state to the left
        end

        if i == Ncells %if last cell
            u_r = u_j; %use same state
        else
            u_r = u(i+1,:); %use state to the left
        end
        
       [F_F(i,:)] = russonov(u_l,u_j,u_r, v_n, mi);
       %[F_F(i,:), ~, ~, ~] = HLLE_flux(u_l,u_j,u_r, v_n, mi, i);
       
       [S(i,:), S_area(i,:), P(i)] = source(u_j, r(i), A(i), dAdx(i), dx, u, Br, P_e_0, i);
         
    end
      
    R = dt'.*((F_F./dx)+S_area+S);
    
    %update state
    u_next = u + R;  

    %% Check for convergence
    SumR = 0;
    nR = size(R,1)*4; %total number of entries in R
    R_list = reshape(R./u_next,[nR,1]);
   
    for i = 1:nR
        SumR = SumR + abs(R_list(i))^2; %Cumulative sum of the elements in R_list squared
    end
    
    Rnorm(N) = sqrt(SumR);
   
    %check for error
    if any( u_next(:,[1 2 4]) < 0 ,'all') %throw error if non physical state present
       dt = dt*0.1;
       continue
    else %no error
        u = u_next; %full state update
    end

    %% Check if state output is due
    if count >= count_max
        % Write out restart file
       
        % outputstate = fopen('Z.txt','w');
        % fprintf(outputstate,'%23.16E %23.16E %23.16E %23.16E %23.16E\n',[x u]');
        % fclose(outputstate);
        
        %reset counter
        count = 0;
        %let user now code is still processing
        fprintf('Working... (N = %d, Rnorm = %.2e) \n',N,Rnorm(N));

        %dt = dt_max;
    end

end

%% Write out final restart file

% outputstate = fopen('U.txt','w');
% fprintf(outputstate,'%23.16E %23.16E %23.16E %23.16E %23.16E\n',[x u]');
% fclose(outputstate);

figure(1); clf; hold on
tiledlayout(3,2);
nexttile([1 2]);
    semilogy(Rnorm,'k','Linewidth',1); ylabel('Residual Norm'); xlabel('iteration')
nexttile();
    plot(x_cell,u(:,1),'k','Linewidth',1); ylabel('neutral density')
nexttile();
    plot(x_cell,u(:,2),'k','Linewidth',1); ylabel('plasma density')
nexttile();
    plot(x_cell,u(:,3)./u(:,2),'k','Linewidth',1); ylabel('Ion speed')
nexttile();   
    pe = u(:,4)/1.5;
    plot(x_cell, pe./(e*u(:,2)),'k','Linewidth',1); ylabel('Electron Temperature')

%% Calculate Performance Metrics

m_dot_i = mi*u(end,3)*pi*r(end)^2;
m_dot = mi*u(1,3)*pi*r_0^2 + mi*u(1,1)*v_n*pi*r_0^2;
Power = trapz(x,P);
eta_m = m_dot_i./m_dot;
Thrust = m_dot_i*u(end,3)./u(end,2);
Isp = u(end,3)/(9.81*u(end,2));
eta = Thrust^2/(2*m_dot*Power);

ne = 1/l*trapz(x,u(:,2));
nn = 1/l*trapz(x,u(:,1));
Te = 1/l*trapz(x,u(:,4)./(e*u(:,2)));

%Calculate minimum RMF Current
sigv_i = 2.9e-12* 10 ./ Te.^(3/2);
sigv_n = 6.6e-19 * ((Te/4 - 0.1)./(1 + (Te/4).^1.6)) .* sqrt((8*e*Te)/(pi*me));
nu = ne.*sigv_i + nn.*sigv_n;

A = 100; %factor for << sign conversion to =

%Modified Hugrass penetration condition
Iw = sqrt(A) * (5/4)^(3/2) * r_bar^2 * sqrt((2*pi*f_rmf.*ne.*me.*nu)/(2*mu0));

%Modify Power
Power = Power + Rppu*Iw.^2;
eta = Thrust.^2./(2.*m_dot.*Power);

%Print out performance metrics
fprintf(['\neta = %.2f%%   \n' ...
            'Thrust = %.2f mN \n' ...
            'Isp = %.0f s     \n' ...
            'Power = %.2f kW \n\n'], eta*100, Thrust*1e3, Isp, Power*1e-3);
end

%BEGIN SUBFUNCTIONS -------------------------------------------------------

function [S,S_area, P] = source(u, r, A, dAdx, dx, u_total, Br, p_e_0, i)
%source calculates the source components via the governing equations from 
%the state u. 
%   u = [n_n; n_e; n_e*v_e; 1.5*p_e];
%   F = [n_n*v_n; n_e*v_e; n_e*v_e^2 + P_e/m_i; 5/2*p_e*v_e];

global mi f_c w v_n

%Constants
e = 1.6022e-19; %c
me = 9.109383632e-31; %kg

%Derived variables
v_t = 2/3*w*r; %radius weighted theta electron velocity, 
               % (2*pi/A)*int_0^r (w*r') r' dr' = 2/3 * w * r

%Derived state quantities
pe = u(4)/1.5;
T_eV = pe./(e*u(2)); %eV
T_e = pe./u(2); %J 

%INCLUDE V_THETA FOR THE IONIZATON REACTION RATE

% Ionization Reaction Rate [m^3/s]
avg_speed = sqrt((8*e*T_eV)/(pi()*me));
if T_eV < 5 %eV
    sigv_iz = (1e-20)*((3.97 + 0.643*T_eV - 0.0368*T_eV*T_eV)*(exp(-12.127/T_eV)))*(avg_speed);
else %>= 5 eV
    sigv_iz = (1e-20)*((-1.031e-4*T_eV*T_eV)+(6.386*exp(-12.127/T_eV)))*(avg_speed);
end

%excitation reation rate [m^3/s]
sigv_rad = 1.93e-19 * (exp(-11.6/T_eV))/(sqrt(T_eV)) * avg_speed;

%collision frequencies
ln_lambda = 23 - (0.5*log((u(2)*1e-6)/(T_eV^3)));
nu_ei = (u(2) * 2.9e-12* ln_lambda) ./ (T_eV.^(3/2));

cross_en =  6.6e-19 * ((T_eV/4 - 0.1)./(1 + (T_eV/4).^1.6));
%nu_en = u(1) *sqrt(e*T_eV/me)*cross_en;
nu_en = u(1) *sqrt((8*e*T_eV)/(pi()*me))*cross_en;

%the velocity is different in different books

%inelastic collision energy losses 
   eps_iz = e*12.1; %J
 eps_nrad = e*8.32; %J
 eps_irad = e*14.8; %J
eps_wall = T_e.*(5/2 + 2* log(sqrt((2*mi)/(pi*me))));
 
 % floating_pot = (T_e/e)*log(sqrt((2*mi)/(pi))); %V
 % floating_pot = floating_pot*e; %eV
 % eps_wall = (2*T_eV + floating_pot)*e; %J (pg477)

%reaction rates [1/m^3 1/s]
  R_iz = u(1).*u(2).*sigv_iz;
R_nrad = u(1).*u(2).*sigv_rad;
R_irad = u(2).^2.*sigv_rad;
R_wall = 0.6*u(2).*sqrt(pe./(u(2)*mi))*2/r*f_c;

%Electric fields
v_e = u(3)./u(2);

%E_t = + me/e*v_t*(nu_ei + nu_en); %sign is irrelevant here. 
E_t = me/e*v_t*(nu_en) +v_e*Br; %nu_ei

F = e*u(2)*E_t; % N/m^3
P = 3/4*pi*F*w*r^3; % W/m % factor 3/4 = 9/8 * 2/3 to correct for improper area integral 

S = [ - R_iz + R_wall, ...
     + R_iz - R_wall, ...
     + e*u(2)/mi*v_t*Br, ...
     - R_iz*eps_iz - R_wall*eps_wall - R_nrad*eps_nrad  + 9/8*e*(u(2)*v_t*E_t)]; %+ u(2)*0.5*r^2*nu_ei*me*w^2]; %- R_irad*eps_irad
                                                                           %9/8 to correct for not perfroming area average with v_t^2.
%note power above = joule heating * A;

% S_iz = [ - R_iz, R_iz , 0, -R_iz*eps_iz];
% 
% S_wall = [R_wall, - R_wall, 0, -R_wall*eps_wall]; %also sketchy - R_wall*eps_wall
% 
% S_rad = [0, 0, 0, - R_nrad*eps_nrad - R_irad*eps_irad];
% 
% S_mom = [0,0,e*u(2)/mi*v_t*Br, 0];
% 
% S_ene = [0,0,0,  (9/8)*u(2)*v_t*me*v_t*(nu_ei+nu_en)]; %nu_en is the problem child
% 
% %S =  S_mom +S_ene+S_iz+S_wall+S_rad;
% S = S_iz;

[S_area] = source_area(u,A,dAdx, dx, v_n, u_total, p_e_0, i);

end

%--------------------------------------------------------------------------
function [S_area] = source_area(u, A, dAdx, dx, v_n, u_total, p_e_0, i)
%computes additional "source" due to unequal areas on right and left of
%cells

A_var = 1/A*dAdx;
pe = u(4)/1.5;
pe_total = [p_e_0*1.5; u_total(:,4); u_total(end,4)];
pe_total = pe_total./1.5;
dpe_dx = (pe_total(i+2) - pe_total(i)) ./ (2*dx);

S_area = [-A_var*u(1)*v_n, ...
          -A_var*u(3), ...
          -A_var*u(3)^2./u(2), ... 
          -A_var*(5/2)*pe*u(3)./u(2) + (u(3)./u(2)).*dpe_dx];
end
%END OF CODE --------------------------------------------------------------
    
    
