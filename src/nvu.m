function [f_rhs, u0, misc_ast, T, U] = nvu()
[f_ast, u0_ast, idx_ast, p_ast, misc_ast] = astrocyte_model('P_L', 0.0804);
[f_smc_ec, u0_smc_ec, idx_smc_ec, p_smc_ec] = smc_ec_model('J_PLC', 0.4);
[f_mech, u0_mech, idx_mech, p_mech] = mechanical_model();

na = length(u0_ast);
nse = length(u0_smc_ec);
nmech = length(u0_mech);

i_a = 1:na;
i_se = (1:nse) + na;
i_mech = (1:nmech) + na + nse; 

u0 = [u0_ast; u0_smc_ec; u0_mech];
rhs(0, u0);
f_rhs = @rhs;
tic
odeopts = odeset('MaxStep', 1, 'RelTol', 1e-3, 'AbsTol', 1e-3);
[T, U] = ode15s(f_rhs, [0 500], u0, odeopts);
toc
    function du = rhs(t, u)
        du = zeros(size(u));
        
        % Get state variables for each sub-model
        u_a = u(i_a);
        u_se = u(i_se);
        u_mech = u(i_mech);
        
        % Evaluate mechanical model first
        Ca_i = u_se(idx_smc_ec.Ca_i);
        [du(i_mech), R, h] = f_mech(t, u_mech, Ca_i);

        % Now SMC/EC 
        K_p = u_a(idx_ast.K_p, :);
        [du(i_se), J_KIR_i] = f_smc_ec(t, u_se, R, h, K_p);
        
        % Now the astrocyte model
        du(i_a) = f_ast(t, u_a, J_KIR_i);
    end

end