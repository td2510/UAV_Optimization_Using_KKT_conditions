function [q,tau,eta,g_x,P_un,P_sn] = Linear_model2(tau_j,eta_j,q_j_1,N)

global P_s V_max sigma_sq H delta_t omega_0 P_c alpha miu q_I1 q_F1 w_s ....
    w_d epsilon sigma Euler eta_max S E_tot Theta Theta_0 P_u P_h

Max_Iteration = 3; % maximum number of iteration

P_u1 = P_u.*ones(1,N);
P_s1 = P_s.*ones(1,N);
P_c1 = P_c.*ones(1,N);

iter = 1; err = 1; 
object_array = [];
err_array = [];
q_array = [];
tau_array = [];
eta_array = [];
P_u_array = [];
P_s_array = [];
while ( (iter < Max_Iteration)&&(err > epsilon) )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DTS optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculating E_fly
    delta = 0.012; rho = 1.225; A = 0.8; s = 0.05; Omega = 100; R = 0.08;
    W = 0.5 ; k = 0.1; d_0 = 0.0151/s/A;

    P_0 = delta*rho*s*A*(Omega*R)^3/8; 
    P_1 = (1+k)*W^1.5/sqrt(2*rho*A);
    P_p = 0.5*d_0*rho*s*A; B = 3/(Omega*R)^2; 
    v_0 = sqrt(W/(2*rho*A)); C = 1/(4*v_0^4);  D = sqrt(C);

    q_j_11 = q_j_1(:,[1:N]); 
    q_j_12 = q_j_1(:,[2:N+1]);
    E_fly = zeros(1,N);
    if iter==1
        for n1=1:N
            E_fly(1,n1) = P_0.*(delta_t + B.*sum((q_j_1(:,n1+1)-q_j_1(:,n1)).^2)./delta_t ) + ...
            P_1.*sqrt( (delta_t.^4 + D.^2.*(sum((q_j_1(:,n1+1) - q_j_1(:,n1)).^2)).^2 ).^0.5- D*sum((q_j_1(:,n1+1) - q_j_1(:,n1)).^2) )...
            + P_p.* pow_pos(norm(q_j_1(:,n1+1)-q_j_1(:,n1)),1.5)./(delta_t.^2);
        end
    else
        for n1=1:N
            E_fly(1,n1) = P_0.*(delta_t + B.*sum((q_j_1(:,n1+1)-q_j_1(:,n1)).^2)./delta_t ) + ...
            P_1.*y_n(:,n1)+ P_p.* pow_pos(norm(q_j_1(:,n1+1)-q_j_1(:,n1)),1.5)./(delta_t.^2);
        end
    end

    %% Calculate Rate from source to UAV and from UAV to destination
    %% Since the R_d usually less than R_u due to it only reflects a part of power thus reflecting rate should less than information rate
    R_u = log2(1+ Theta_0.*P_s1./( H^2+sum( (q_j_12 - w_s).^2 )).^(alpha/2) );
%     R_d = log2(1+ Theta.*(eta_j.*omega_0.*P_s1+P_u1.*(ceil(sigma)).*( H^2+sum( (q_j_12 - w_s).^2 )).^(alpha/2))...
%     ./( H^2+sum( (q_j_12 - w_s).^2 )).^(alpha/2) ./( H^2+sum( (q_j_12 - w_d).^2 )).^(alpha/2) );
    R_d = log2(1+ Theta.*(eta_j.*omega_0.*P_s1+P_u1.*(1+ceil(sigma)).*( H^2+sum( (q_j_12 - w_s).^2 )).^(alpha/2))...
    ./( H^2+sum( (q_j_12 - w_s).^2 )).^(alpha/2) ./( H^2+sum( (q_j_12 - w_d).^2 )).^(alpha/2) );

    indice =  find(R_d > R_u); % It should be emply

    d_ns =  ( H^2+sum( (q_j_12 - w_s).^2 )).^(alpha/2); % Distance from UAV to source at time slot n
    Xi_1 = miu.*delta_t.*omega_0.*P_h./d_ns;

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
    %% Updating tau_j
    tau_j = tau;
    tau_j(tau_j >= 1)=0.9;
    tau_j(tau_j <= 0)=0.1;  % In this case the E_fly is larger than EH in that time slot, thus all time slot should used for EH
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Updating eta_j
    eta_j = eta_max; % we dont need to optimize \eta since it is a linear function
    eta = eta_max;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% UAV trajectory optimization for Linear Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z_1j = ( H^2+sum( (q_j_12 - w_s).^2 )).^(alpha/2);
    z_2j = ( H^2+sum( (q_j_12 - w_d).^2 )).^(alpha/2);
    [q,z_1,z_2,y_n,g_x] = Trajectory_Linear2(q_j_11,q_j_12,z_1j,z_2j,tau_j,eta_j,N,P_s1,P_u1);

    q_j_1 = q; % Updateing q_j_1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Checking the convergence condition
    if ~isnan(g_x)
        object_array = [object_array;g_x];
        q_array = [q_array;q];
        tau_array = [tau_array;tau];
        eta_array = [eta_array;eta];
        P_u_array = [P_u_array;P_u1];
        P_s_array = [P_s_array;P_s1];
    end
    if iter>1 && ~isnan(g_x)
        err =  abs(object_array(iter)-object_array(iter-1));
        err_array= [err_array; err];
    end
    iter = iter+1;    
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Update maximum value     
if ~isempty (object_array) % For 3D-2UAV case
    [value,idx] = max(object_array);
    g_x = object_array(idx);
    q = q_array([2*idx-1:2*idx],:);
    tau = tau_array(idx,:);
    eta = eta_array(idx,:);
    P_un = P_u_array(idx,:);
    P_sn = P_s_array(idx,:);
else
    q = NaN; tau = NaN; eta = NaN;
    P_un = NaN; P_sn = NaN;
end  

end