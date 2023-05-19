function [q,z_1,z_2,y_n,g_x] = Trajectory_Linear_2(q_j_21,q_j_22,z_1j,z_2j,tau_j,eta_j,N,P_s1,P_u1)

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
P_1 = (1+k)*W^1.5/sqrt(2*rho*A);
P_p = 0.5*d_0*rho*s*A; B = 3/(Omega*R)^2; 
v_0 = sqrt(W/(2*rho*A)); C = 1/(4*v_0^4);  D = sqrt(C);
P_c1 = P_c.*ones(1,N);

y_j = sqrt( (delta_t.^4 + D.^2.*(sum((q_j_22 - q_j_21).^2)).^2 ).^0.5- D*sum((q_j_22 - q_j_21).^2) );%44
while ((err>epsilon)&&(iter<=3))

%% Run the CVX to solve the problem P3.2
cvx_begin %quiet
   variable q(3,N+1) %???
   variable z_1(1,N)
   variable z_2(1,N)
   variable y_n(1,N)
%    variable E_fly(1,N)
   
   w_s1 = w_s.*ones(3,N);%???
   w_d1 = w_d.*ones(3,N);%???
   z1_Taylor = 1./z_1j - (z_1-z_1j)./square(z_1j);
   
   bar_P_u1 = P_u1*(1+ceil(sigma));
   
   Theta_1 = log2(1+Theta_0.*P_s1./z_1j)-Theta_0.*P_s1.*(z_1-z_1j)./z_1j./log(2)./(z_1j+Theta_0.*P_s1);
   Theta_2 = log2(1+Theta.*(eta_j.*omega_0.*P_s1 + bar_P_u1.*z_1j)./(z_1j.*z_2j) )...
        -Theta.*eta_j.*omega_0.*P_s1.*(z_1-z_1j)./z_1j./log(2)./(Theta.*eta_j.*omega_0.*P_s1+z_1j.*(z_2j+Theta.*bar_P_u1))...
        -Theta.*(eta_j.*omega_0.*P_s1+bar_P_u1.*z_1j).*(z_2-z_2j)./z_2j./log(2)./(Theta.*eta_j.*omega_0.*P_s1+z_1j.*(z_2j+Theta.*bar_P_u1));
%    Theta_3 =  log2(1+Theta_0.*P_u1./z_2j)-Theta_0.*P_u1.*(z_2-z_2j)./z_2j./log(2)./(z_2j+Theta_0.*P_u1); % Data transmission rate from UAV to d if it cached a part of f file 
   g_x = sum(tau_j.*delta_t.*Theta_2 ); %(49a)
%       g_x = sum(tau_j.*delta_t.*Theta_2 + ceil(sigma).*tau_j.*delta_t.*Theta_3 );
%     if sigma == 0
%         g_x = sum(tau_j.*delta_t.*Theta_2 );
%     else
%         g_x = sum(tau_j.*delta_t.*Theta_2 + tau_j.*delta_t.*Theta_3 );
%     end
   %% Maximization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   minimize (-g_x)  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   subject to   
    ii = 1: N;
    norm(q(:,ii+1)-q(:,ii)) <= V_max.*delta_t; % constraint 22e
    q(:,1)==q_I2; q(:,N+1)==q_F2; % Constraint 22f        
    sum( (q(:,ii+1)-w_d1).^2 ) <= (z_2).^(2/alpha); % Constraint 37c %???
    sum( (q(:,ii+1)-w_s1).^2 ) <= (z_1).^(2/alpha); % Constraint 37b, sum( (q(:,ii+1)-w_s1).^2 ) %???
    sum(tau_j.*delta_t.*Theta_1) + sigma.*S >= g_x; % Constraint 43c
%     sum(tau_j.*delta_t.*Theta_1) >= S.*(1-sigma); % Constraint 53c
    g_x >= S; % Constraint 43d
    q(2,:) >= 0; % Keep for y value >= 0
    q(3,:) >= 3; %???
    q(3,:) <= 10; %???
    %% energy constraint
    E_fly = cvx(zeros(1,N));
    for n=1:N
        delta_t^4.*pow_p(y_n(:,n),-2)<=  pow_p(y_j(:,n),2)+ 2.*y_j(:,n).*(y_n(:,n)-y_j(:,n))...
        -2.*D.*sum( (q_j_22(:,n)-q_j_21(:,n)).^2 )+ 4.*D.*(q_j_22(:,n)-q_j_21(:,n))'*(q(:,n+1)-q(:,n));%49c
        E_fly(1,n) = P_0.*(delta_t + B.*sum((q(:,n+1)-q(:,n)).^2) ) + ...
        P_1.*y_n(:,n)+P_p.* pow_pos(norm(q(:,n+1)-q(:,n)),1.5)./(delta_t.^2);
    end
    
    for n=1:N               
        tau_j1 = tau_j(:,[1:n]);
        Sum_E_Fly = E_fly(:,[1:n]);
        z1_Taylor_n = z1_Taylor(:,[1:n]);
        P_c1n = P_c1(:,[1:n]); %???
        P_u1n = P_u1(:,[1:n]); %???
        sum(Sum_E_Fly+tau_j1.*delta_t.*(P_c1n+P_u1n))<=sum(miu.*(1-tau_j1).*delta_t.*omega_0.*P_h.*z1_Taylor_n);% Constraint 47c
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

if ~isnan(q)   
    time = toc(start_CVX);
else
    q = NaN;
    time = inf;
end

end