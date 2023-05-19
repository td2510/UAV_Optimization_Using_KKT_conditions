function [q,tau,eta,g_x] = Linear_model_1(tau_j,eta_j,q_j_1,N)

global P_s V_max sigma_sq H delta_t omega_0 P_c alpha miu q_I1 q_F1 w_s ....
    w_d epsilon sigma Euler eta_max S E_tot Theta Theta_0 P_u P_h

iter = 1; err = 1; iter2 = 1;
Max_Iteration = 5; % maximum number of iteration
object_array = [];
err_array = [];

q_array = [];
tau_array = [];
eta_array = [];

P_u1 = P_u.*ones(1,N);
P_s1 = P_s.*ones(1,N);
P_c1 = P_c.*ones(1,N);
while( iter2 < 4)
    iter = 1; err = 1; 
while ( (iter < Max_Iteration)&&(err > epsilon) )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DTS optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Calculating E_fly
    delta = 0.012; rho = 1.225; A = 0.8; s = 0.05; Omega = 100; R = 0.08;
    W = 0.5 ; k = 0.1; d_0 = 0.0151/s/A;

    P_0 = delta*rho*s*A*(Omega*R)^3/8; 
    P_1 = (1+k)*W^1.5/sqrt(2*rho*A); %k = I
    P_p = 0.5*d_0*rho*s*A; B = 3/(Omega*R)^2; % B = k1, P_p = k3
    v_0 = sqrt(W/(2*rho*A)); C = 1/(4*v_0^4);  D = sqrt(C); % C =  k2^2
    
    q_j_11 = q_j_1(:,[1:N]);  %q_j_11 = 1:(length(q_j)-1);
    q_j_12 = q_j_1(:,[2:N+1]);%q_j_12 = 2:length(q_j);
    if iter==1
        E_fly = P_0.*(delta_t + B.*sum((q_j_12 - q_j_11).^2) ) + ...
        P_1.*sqrt( (delta_t.^4 + D.^2.*sum((q_j_12 - q_j_11).^2).^2 ).^0.5- D*sum((q_j_12 - q_j_11).^2) )...
        + P_p.* pow_pos(norm(q_j_12-q_j_11),1.5)./(delta_t.^2)
    else
        E_fly = P_0.*(delta_t + B.*sum((q_j_12 - q_j_11).^2) ) + ...
        P_1.*y_n + P_p.* pow_pos(norm(q_j_12-q_j_11),1.5)./(delta_t.^2);
        
    end
   
    %% Calculate Rate from source to UAV and from UAV to destination
    %% Since the R_d usually less than R_u due to it only reflects a part of power thus reflecting rate should less than information rate
    R_u = log2(1+ Theta_0.*P_s1./( sum( (q_j_12 - w_s).^2 )).^(alpha/2) ); %???
    R_d = log2(1+ Theta.*(eta_j.*omega_0.*P_s1+P_u1*(1+ceil(sigma)).*( sum( (q_j_12 - w_s).^2 )).^(alpha/2))./( sum( (q_j_12 - w_s).^2 )).^(alpha/2) ./( sum( (q_j_12 - w_d).^2 )).^(alpha/2) );%???
    %P_u = P_u*(1+ceil(sigma))
    indice =  find(R_d > R_u); % It should be empty
    
    d_ns =  ( sum( (q_j_12 - w_s).^2 )).^(alpha/2); % Distance from UAV to source at time slot n %???
    Xi_1 = miu.*delta_t.*omega_0.*P_h./d_ns

    tau = zeros(1, N);
    if ~isempty(indice)
        for i = 1: N
            if ~isempty(find(indice == i))
                tau(1,i)= sigma*S/N/(R_d(i)-R_u(i));
            else
                tau(1,i)= (Xi_1(1,i) - E_fly(1,i))./(Xi_1(1,i)+delta_t.*(P_c1(1,i)+P_u1(1,i)));
            end
        end
    else
        tau = (Xi_1 - E_fly)./(Xi_1+delta_t.*(P_c1+P_u1));
    end 
    tau(tau>1)=0.9;
    tau(tau<=0)=0.1;
    %% Updating tau_j
    tau_j = tau;
    idx = find(tau<0); % In this case the E_fly is larger than EH in that time slot, thus all time slot should used for EH
    tau_j(idx)=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Updating eta_j
    eta_j = eta_max; % we dont need to optimize \eta since it is a linear function
    eta = eta_max;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% UAV trajectory optimization for Linear Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z_1j = ( sum( (q_j_12 - w_s).^2 )).^(alpha/2); %???
    z_2j = ( sum( (q_j_12 - w_d).^2 )).^(alpha/2); %???
    
    [q,z_1,z_2,y_n,g_x] = Trajectory_Linear_1(q_j_11,q_j_12,z_1j,z_2j,tau_j,eta_j,N,P_s1,P_u1);
    
    q_j_1 = q; % Updateing q_j
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Checking the convergence condition
    bar_P_u = P_u1*(1+ceil(sigma));
    Theta_2 = log2(1+Theta.*(eta_j.*omega_0.*P_s1 + bar_P_u.*z_1j)./(z_1j.*z_2j) )...
        -Theta.*eta_j.*omega_0.*P_s1.*(z_1-z_1j)./z_1j./log(2)./(Theta.*eta_j.*omega_0.*P_s1+z_1j.*(z_2j+Theta.*bar_P_u))...
        -Theta.*(eta_j.*omega_0.*P_s1+bar_P_u.*z_1j).*(z_2-z_2j)./z_2j./log(2)./(Theta.*eta_j.*omega_0.*P_s1+z_1j.*(z_2j+Theta.*bar_P_u));
    g_x = sum(tau_j.*delta_t.*Theta_2);
    
    if ~isnan(g_x)
        object_array = [object_array;g_x];
        q_array = [q_array;q];
        tau_array = [tau_array;tau];
        eta_array = [eta_array;eta];
    end
    if iter>1 && ~isnan(g_x)
        err =  abs(object_array(iter)-object_array(iter-1));
        err_array= [err_array; err];
    end
    iter = iter+1;    
end
%% Power optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [P_un,P_sn,g_x] = Optimize_P_1(P_u1,P_s1,tau_j,eta_j,N,q_j_1);
    P_u1 = P_un;
    P_s1 = P_sn;
iter2 = iter2+1;
end
%     if ~isempty (object_array)
%         [value,idx] = max(object_array);
%         g_x = object_array(idx);
%         q = q_array([2*idx-1:2*idx],:);
%         tau = tau_array(idx,:);
%         eta = eta_array(idx,:);
%     else
%         q = NaN; tau = NaN; eta = NaN;
%     end

end