function [P_un,P_sn,g_x]= Optimize_P_2(P_u1,P_s1,tau_j,eta_j,N,q_j_2)

global P_s V_max sigma_sq H delta_t omega_0 P_c alpha miu q_I2 q_F2 w_s ....
    w_d epsilon sigma Euler eta_max S E_tot Theta Theta_0 P_u P_h
start_CVX = tic;

err = 1.0;

iter = 0;
err_matrix = []; % Including all the err values
gx_array = []; % All objective values

% cvx_solver mosek

delta = 0.012; rho = 1.225; A = 0.8; s = 0.05; Omega = 100; R = 0.08;
W = 0.5 ; k = 0.1; d_0 = 0.0151/s/A;

P_0 = delta*rho*s*A*(Omega*R)^3/8; 
P_1 = (1+k)*W^1.5/sqrt(2*rho*A); %k = I
P_p = 0.5*d_0*rho*s*A; B = 3/(Omega*R)^2; % B = k1, P_p = k3
v_0 = sqrt(W/(2*rho*A)); C = 1/(4*v_0^4);  D = sqrt(C); % C =  k2^2
q_j_21 = q_j_2(:,[1:N]);
q_j_22 = q_j_2(:,[2:N+1]);


y_j = sqrt( (delta_t.^4 + D.^2.*(sum((q_j_22 - q_j_21).^2)).^2 ).^0.5- D*sum((q_j_22 - q_j_21).^2) );%44
while ((err>epsilon)&&(iter<=5))
%% Run the CVX to solve the problem P3.2
cvx_begin %quiet
   variable P_un(1,N)
   variable P_sn(1,N)
   
   P_c1 = P_c.*ones(1,N);

   bar_P_u = P_u1*(1+ceil(sigma));
   z_1j = ( sum( (q_j_22 - w_s).^2 )).^(alpha/2); %???
   z_2j = ( sum( (q_j_22 - w_d).^2 )).^(alpha/2); %???
    
   Theta_1 = log2(1+Theta_0.*P_s1./z_1j)+...
       (Theta_0.*(P_sn-P_s1)./log(2))./(z_1j+Theta_0.*P_s1);
%    Theta_2 = log2(1+Theta.*(eta_j.*omega_0.*P_s1 + bar_P_u.*z_1j)./(z_1j.*z_2j) )...
%         +(Theta.*eta_j.*omega_0.*(P_sn-P_s1)+Theta.*z_1j.*(ceil(sigma)).*(P_un-P_u1))...
%         ./log(2)./((z_1j.*z_2j)+Theta.*eta_j.*omega_0.*P_s1+Theta.*z_1j.*bar_P_u);
   Theta_2 = log2(1+Theta.*(eta_j.*omega_0.*P_s1 + bar_P_u.*z_1j)./(z_1j.*z_2j) )...
        +(Theta.*eta_j.*omega_0.*(P_sn-P_s1)+Theta.*z_1j.*(1+ceil(sigma)).*(P_un-P_u1))...
        ./log(2)./((z_1j.*z_2j)+Theta.*eta_j.*omega_0.*P_s1+Theta.*z_1j.*bar_P_u);
   g_x = sum(tau_j.*delta_t.*Theta_2 ); %(49a)
%    E_fly = P_0.*(delta_t + B.*sum((q_j_22 - q_j_21).^2) ) + ...
%         P_1.*y_j+P_p.* pow_pos(norm(q_j_22 - q_j_21),1.5)./(delta_t.^2);
   E_fly = zeros(1,N);
    for n1=1:N
        E_fly(1,n1) = P_0.*(delta_t + B.*sum((q_j_2(:,n1+1)-q_j_2(:,n1)).^2) ) + ...
        P_1.*y_j(:,n1)+P_p.* pow_pos(norm(q_j_2(:,n1+1)-q_j_2(:,n1)),1.5)./(delta_t.^2);
    end
   %% Maximization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   minimize (-g_x)  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   subject to   
    sum(tau_j.*delta_t.*Theta_1) + sigma.*S >= g_x; 
    g_x >= S; 
    P_un >= 0;
    P_un <= 4;
    P_sn >= 0;
    (P_un+P_sn) <= 100;
    (P_un+P_sn) >= 0;
    %% energy constraint
%      sum(E_fly+tau_j.*delta_t.*(P_c1+P_un)) <= sum((miu.*(1-tau_j).*delta_t.*omega_0.*P_h)./z_1j);
    for n=1:N               
        tau_j1 = tau_j(:,[1:n]);
        Sum_E_Fly = E_fly(:,[1:n]);
        z_1jn = z_1j(:,[1:n]);
        P_c1n = P_c1(:,[1:n]); %???
        P_unn = P_un(:,[1:n]); %???
        sum(Sum_E_Fly+tau_j1.*delta_t.*(P_c1n+P_unn))<=sum((miu.*(1-tau_j1).*delta_t.*omega_0.*P_h))./z_1jn;% Constraint 47c
    end
    %%
cvx_end
iter = iter + 1;
gx_array = [gx_array,g_x];
if iter>1 && ~isnan(g_x)
    err =  abs(gx_array(iter)-gx_array(iter-1));
    err_matrix = [err_matrix; err];
end

end % End while
%% Calculate the total energy consumtion

if ~isnan(P_un)   
    time = toc(start_CVX);
else
    P_un = NaN;
    time = inf;
end
if ~isnan(P_sn)   
    time = toc(start_CVX);
else
    P_sn = NaN;
    time = inf;
end

end