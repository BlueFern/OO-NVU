classdef ANLS < handle

    properties
        params
        u0
        index
        n_out
        idx_out
        enabled
    end
    methods
        function self = ANLS(varargin)
            self.params = parse_inputs(varargin{:});
            self.index = indices();
            self.u0 = initial_conditions(self.index);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices();
        end
        function [du, varargout] = rhs(self, t, u, R, J_pump1_sa, K_e, Glu, K_s, Na_e)
            % Initalise inputs and parameters
            idx = self.index;
            p = self.params;
            t = t(:).';
            
            GLUe = (Glu ./ 1000); %uM to mM
            
            GLCn = u(idx.GLCn, :);            
            G6Pn = u(idx.G6Pn, :);
            F6Pn = u(idx.F6Pn, :);
            GAPn = u(idx.GAPn, :);
            PEPn = u(idx.PEPn, :);
            PYRn = u(idx.PYRn, :);
            LACn = u(idx.LACn, :);
            NADHn = u(idx.NADHn, :);
            ATPn = u(idx.ATPn, :);
            GLCg = u(idx.GLCg, :);
            F6Pg = u(idx.F6Pg, :);
            G6Pg = u(idx.G6Pg, :);
            GLYg = u(idx.GLYg, :);  
            GLCe = u(idx.GLCe, :);
            LACe = u(idx.LACe, :);
            O2c = u(idx.O2c, :);
            GLCc = u(idx.GLCc, :);
            LACc = u(idx.LACc, :);
            LACg = u(idx.LACg, :);
            O2n = u(idx.O2n, :);
            GAPg = u(idx.GAPg, :);
            PEPg = u(idx.PEPg, :);
            PYRg = u(idx.PYRg, :);
            NADHg = u(idx.NADHg, :);
            O2g = u(idx.O2g, :);
            PCrg = u(idx.PCrg, :);
            ATPg = u(idx.ATPg, :);
            PCrn = u(idx.PCrn, :);
            
            GLUg = u(idx.GLUg, :);
            Nag = u(idx.Nag, :);
            %GLUn = u(idx.GLUn, :);
            %GLUe = u(idx.GLUe, :);
            
            
            
%             %CALC algebraic vars
%             V_en_GLC =  p.Vm_en_GLC .* (GLCe ./ (GLCe + p.Km_en_GLC) - GLCn ./ (GLCn + p.Km_en_GLC));
%             Vn_hk =  p.Vmax_n_hk .* ATPn .* (GLCn ./ (GLCn + p.Km_GLC)) .* (1.0 - 1.0 ./ (1.0 + exp(  - p.aG6P_inh_hk .* ( 1.0 .* (G6Pn - p.G6P_inh_hk)))));
%             Vn_pgi =  p.Vmaxf_n_pgi .* (G6Pn ./ (G6Pn + p.Km_G6P)) -  p.Vmaxr_n_pgi .* (F6Pn ./ (F6Pn + p.Km_F6P_pgi));
%             Vn_pfk =  p.kn_pfk .* ATPn .* (F6Pn ./ (F6Pn + p.Km_F6P_pfk)) .* power(1.0 + power(ATPn ./ p.Ki_ATP, p.nH), -1.0);
%             Vcg_GLC =  p.Vm_cg_GLC .* (GLCc ./ (GLCc + p.Km_cg_GLC) - GLCg ./ (GLCg + p.Km_cg_GLC));
%             Veg_GLC =  p.KO1 .* p.Vm_eg_GLC .* (GLCe ./ (GLCe + p.Km_eg_GLC) - GLCg ./ (GLCg + p.Km_eg_GLC));
%             Vg_hk =  p.Vmax_g_hk .* ATPg .* (GLCg ./ (GLCg + p.Km_GLC)) .* (1.0 - 1.0 ./ (1.0 + exp(  - p.aG6P_inh_hk .* ( 1.0 .* (G6Pg - p.G6P_inh_hk)))));
%             Vg_pgi =  p.Vmaxf_g_pgi .* (G6Pg ./ (G6Pg + p.Km_G6P)) -  p.Vmaxr_g_pgi .* (F6Pg ./ (F6Pg + p.Km_F6P_pgi));
%             Vg_pfk =  p.kg_pfk .* ATPg .* (F6Pg ./ (F6Pg + p.Km_F6P_pfk)) .* power(1.0 + power(ATPg ./ p.Ki_ATP, p.nH), -1.0);
%             Vg_glys =  p.Vmax_glys .* (G6Pg ./ (G6Pg + p.Km_G6P_glys)) .* (1.0 - 1.0 ./ (1.0 + exp(  - p.aGLY_inh .* ( 1.0 .* (GLYg - p.GLY_inh)))));
%             unitstepSB2 = piecewise({VOI - (p.tend_GLY + p.to + p.to_GLY)>=0.0, 1.0 }, 0.0);
%             deltaVt_GLY = 1.0 +  p.stim .* ( p.delta_GLY .* p.KO3 .* (1.0 ./ (1.0 + exp( 1.0 .*  - p.sr_GLY .* (VOI - (p.to + p.to_GLY))))) .* (1.0 - unitstepSB2));
%             Vg_glyp =  p.Vmax_glyp .* (GLYg ./ (GLYg + p.Km_GLY)) .* deltaVt_GLY;
%             Vce_GLC =  p.Vm_ce_GLC .* (GLCc ./ (GLCc + p.Km_ce_GLC) - GLCe ./ (GLCe + p.Km_ce_GLC));
%             Vne_LAC =  p.Vmax_ne_LAC .* (LACn ./ (LACn + p.Km_ne_LAC) - LACe ./ (LACe + p.Km_ne_LAC));
%             Vge_LAC =  p.Vmax_ge_LAC .* (LACg ./ (LACg + p.Km_ge_LAC) - LACe ./ (LACe + p.Km_ge_LAC));
%             Vec_LAC =  p.Vm_ec_LAC .* (LACe ./ (LACe + p.Km_ec_LAC) - LACc ./ (LACc + p.Km_ec_LAC));
%             Vg_leak_Na =  (p.Sm_g ./ p.Vg) .* (p.gg_NA ./ F) .* ( (RT ./ F) .* log(NAe ./ NAg) - Vm);
%             Vg_pump =  (p.Sm_g ./ p.Vg) .* kpump .* ATPg .* NAg .* power(1.0 + ATPg ./ Km_pump, -1.0);
%             Veg_GLU =  p.Vmax_eg_GLU .* (GLUe ./ (GLUe + p.Km_GLU));
%             Vg_gs =  p.Vmax_g_gs .* ( (GLUg ./ (GLUg + p.Km_GLU)) .* (ATPg ./ (ATPg + p.Km_ATP)));
%             Vcn_O2 =  (p.PScapn ./ p.Vn) .* ( p.Ko2 .* power(p.HbOP ./ O2c - 1.0, -1.0 ./ p.nh_O2) - O2n);
%             Vcg_O2 =  (p.PScapg ./ p.Vg) .* ( p.Ko2 .* power(p.HbOP ./ O2c - 1.0, -1.0 ./ p.nh_O2) - O2g);
%             Fin_t = CBF0 + ( p.stim .* CBF0 .* p.deltaf .* (1.0 ./ (1.0 + exp( ( 1.0 .*  - p.sr) .* (VOI - ((p.to + p.t1) - 3.0))))) -  p.stim .* CBF0 .* p.deltaf .* (1.0 ./ (1.0 + exp( ( 1.0 .*  - p.sr) .* (VOI - (p.to + tend + p.t1 + 3.0))))));
%             Vc_O2 =  2.0 .* (Fin_t ./ p.Vc) .* (p.O2a - O2c);
%             Vc_GLC =  2.0 .* (Fin_t ./ p.Vc) .* (p.GLCa - GLCc);
%             Vgc_LAC =  p.Vmax_gc_LAC .* (LACg ./ (LACg + p.Km_gc_LAC) - LACc ./ (LACc + p.Km_gc_LAC));
%             Vc_LAC =  2.0 .* (Fin_t ./ p.Vc) .* (p.LACa - LACc);
%             Fout_t =  CBF0 .* ((power(Vv ./ p.Vv0, 2.0) +  p.tv .* power(Vv ./ p.Vv0, -0.50) .* (Fin_t ./ p.Vv0)) ./ (1.0 +  CBF0 .* p.tv .* power(Vv ./ p.Vv0, -0.50) .* (1.0 ./ p.Vv0)));
%             Vn_leak_Na =  (p.Sm_n ./ p.Vn) .* (p.gn_NA ./ F) .* ( (RT ./ F) .* log(NAe ./ NAn) - Vm);
%             Vn_pump =  (p.Sm_n ./ p.Vn) .* kpump .* ATPn .* NAn .* power(1.0 + ATPn ./ Km_pump, -1.0);
%             
%             v_stim =  p.stim .* (p.v1_n +  p.v2_n .* ((VOI - p.to) ./ p.t_n_stim) .* exp( - ( (VOI - p.to) .* (unitpulseSB ./ p.t_n_stim)))) .* unitpulseSB;
%             Vn_stim = v_stim;
%             NADn = p.NADH_n_tot - NADHn;
%             Vn_ldh =  p.kfn_ldh .* PYRn .* NADHn -  p.krn_ldh .* LACn .* NADn;
%             Vn_stim_GLU =  Vn_stim .* p.R_GLU_NA .* p.KO2 .* (GLUn ./ (GLUn + p.Km_GLU));
%             NADg = p.NADH_g_tot - NADHg;
%             Vg_ldh =  p.kfg_ldh .* PYRg .* NADHg -  p.krg_ldh .* LACg .* NADg;
%             ADPn =  (ATPn ./ 2.0) .* ( - p.qak + power((power(p.qak, 2.0) +  4.0 .* p.qak .* (p.ATPtot ./ ATPn - 1.0)), 1.0  ./  2));
%             Vn_pgk =  p.kn_pgk .* GAPn .* ADPn .* (NADn ./ NADHn);
%             Vn_pk =  p.kn_pk .* PEPn .* ADPn;
%             Vn_mito =  p.Vmax_n_mito .* (O2n ./ (O2n + p.Km_O2)) .* (ADPn ./ (ADPn + p.Km_ADP)) .* (PYRn ./ (PYRn + p.Km_PYR)) .* (1.0 - 1.0 ./ (1.0 + exp(  - p.aATP_mito .* ( 1.0 .* (ATPn ./ ADPn -  1.0 .* p.rATP_mito)))));
%             CRn = PCrn_tot - PCrn;
%             Vn_ck =  p.kfn_ck .* PCrn .* ADPn -  p.krn_ck .* CRn .* ATPn;
%             Vn_ATPase =  p.Vmax_n_ATPase .* (ATPn ./ (ATPn + 0.0010));
%             u_n = power(p.qak, 2.0) +  4.0 .* p.qak .* (p.ATPtot ./ ATPn - 1.0);
%             dAMP_dATPn = (p.qak ./ 2.0 +  p.qak .* (p.ATPtot ./ ( ATPn .* power(u_n, 1.0  ./  2)))) - (1.0 +  0.50 .* power(u_n, 1.0  ./  2));
%             ADPg =  (ATPg ./ 2.0) .* ( - p.qak + power((power(p.qak, 2.0) +  4.0 .* p.qak .* (p.ATPtot ./ ATPg - 1.0)), 1.0  ./  2));
%             Vg_pgk =  p.kg_pgk .* GAPg .* ADPg .* (NADg ./ NADHg);
%             Vg_pk =  p.kg_pk .* PEPg .* ADPg;
%             Vg_mito =  p.Vmax_g_mito .* (O2g ./ (O2g + p.Km_O2)) .* (ADPg ./ (ADPg + p.Km_ADP)) .* (PYRg ./ (PYRg + p.Km_PYR)) .* (1.0 - 1.0 ./ (1.0 + exp( 1.0 .*  - p.aATP_mito .* (ATPg ./ ADPg -  1.0 .* p.rATP_mito))));
%             CRg = p.PCrg_tot - PCrg;
%             Vg_ck =  p.kfg_ck .* PCrg .* ADPg -  p.krg_ck .* CRg .* ATPg;
%             Vnc_CO2 =  3.0 .* Vn_mito;
%             Vc_CO2 =  2.0 .* (Fin_t ./ p.Vc) .* (CO2c - p.CO2a);
%             Vgc_CO2 =  3.0 .* Vg_mito;
%             Vg_ATPase =  p.Vmax_g_ATPase .* (ATPg ./ (ATPg + 0.0010));
%             u_g = power(p.qak, 2.0) +  4.0 .* p.qak .* (p.ATPtot ./ ATPg - 1.0);
%             dAMP_dATPg = (p.qak ./ 2.0 +  p.qak .* (p.ATPtot ./ ( ATPg .* power(u_g, 1.0  ./  2)))) - (1.0 +  0.50 .* power(u_g, 1.0  ./  2));
%             BOLD =  p.Vv0 .* ( (p.k1 + p.k2) .* (1.0 - p.dHb ./ p.dHb0) -  (p.k2 + p.k3) .* (1.0 - Vv ./ p.Vv0));
%             unitstepSB = piecewise({VOI - (tend + p.to)>=0.0, 1.0 }, 0.0);
%             AMPn = p.ATPtot - (ATPn + ADPn);
%             AMPg = p.ATPtot - (ATPg + ADPg);
            
            
            % CALC RATES
            CBF2 = p.CBF_init * (R.^4 / p.R_init^4);
            
            pumpScaling = 2.232; %2.25 for low GLC, 2.232
            atp_term = p.ATPheight * 0.5 * (1 + tanh(((ATPn - p.ATPshift) / p.ATPslope))); % new term trying to drop pump
            atp_term2 = ((1 + (p.ATP_init_n ./ ATPn)).^(-1)); %term from chang supp material
            I_pumpv2 = pumpScaling * p.Imax * J_pump1_sa .* atp_term2; 
            I_pump = pumpScaling * p.Imax * J_pump1_sa .* atp_term;
            Vn_pump = 6.6 .* (p.As / p.Vs) * 1e2 * I_pumpv2 ./ p.F;
            %Vn_pump = 0.1583 + (0.035 * rectpuls(t - 102.5, 5));
            
            K_s_mM = K_s ./ 1000; 
            Jk_pump  = (1 + (p.K_init_s ./ K_s_mM)).^(-2) .* (1 + (p.Na_init_k ./ Nag)) .^ (-3) .* ((1 + (p.ATP_init_k ./ ATPg)).^(-1));
            Ik_pump = p.Imax_k .* Jk_pump; 
            Vk_pump = 3 .* (1 / p.R_k) * Ik_pump ./ p.F; %3.7
            
            
            
            Vce_GLC =  p.Vm_ce_GLC .* (GLCc ./ (GLCc + p.Km_ce_GLC) - GLCe ./ (GLCe + p.Km_ce_GLC));
            Vcg_GLC =  p.Vm_cg_GLC .* (GLCc ./ (GLCc + p.Km_cg_GLC) - GLCg ./ (GLCg + p.Km_cg_GLC));
            
                             
            V_en_GLC =  p.Vm_en_GLC .* (GLCe ./ (GLCe + p.Km_en_GLC) - GLCn ./ (GLCn + p.Km_en_GLC));
            Vn_hk =  p.Vmax_n_hk .* ATPn .* (GLCn ./ (GLCn + p.Km_GLC)) .* (1.0 - 1.0 ./ (1.0 + exp(  - p.aG6P_inh_hk .* ( 1.0 .* (G6Pn - p.G6P_inh_hk)))));
            Vn_pgi =  p.Vmaxf_n_pgi .* (G6Pn ./ (G6Pn + p.Km_G6P)) -  p.Vmaxr_n_pgi .* (F6Pn ./ (F6Pn + p.Km_F6P_pgi));
            Vn_pfk =  p.kn_pfk .* ATPn .* (F6Pn ./ (F6Pn + p.Km_F6P_pfk)) .* power(1.0 + power(ATPn ./ p.Ki_ATP, p.nH), -1.0);
            Veg_GLC =  p.KO1 .* p.Vm_eg_GLC .* (GLCe ./ (GLCe + p.Km_eg_GLC) - GLCg ./ (GLCg + p.Km_eg_GLC));
            Vg_hk =  p.Vmax_g_hk .* ATPg .* (GLCg ./ (GLCg + p.Km_GLC)) .* (1.0 - 1.0 ./ (1.0 + exp(  - p.aG6P_inh_hk .* ( 1.0 .* (G6Pg - p.G6P_inh_hk)))));
            Vg_pgi =  p.Vmaxf_g_pgi .* (G6Pg ./ (G6Pg + p.Km_G6P)) -  p.Vmaxr_g_pgi .* (F6Pg ./ (F6Pg + p.Km_F6P_pgi));
            Vg_pfk =  p.kg_pfk .* ATPg .* (F6Pg ./ (F6Pg + p.Km_F6P_pfk)) .* power(1.0 + power(ATPg ./ p.Ki_ATP, p.nH), -1.0);

            Vg_glys =  p.Vmax_glys .* (G6Pg ./ (G6Pg + p.Km_G6P_glys)) .* (1.0 - 1.0 ./ (1.0 + exp(  - p.aGLY_inh .* ( 1.0 .* (GLYg - p.GLY_inh)))));
            
            unitstepSB2 = 0.0;
            % For GLY
            if t - (p.t_0 + p.lengthpulse + p.tend_GLY) >= 0.0
                unitstepSB2 = 1.0;
            end
            
            gly_switch = 0;
            
            if p.lengthpulse > p.tstart_GLY
                gly_switch = 1;
            end
            
            %unitstepSB2 = piecewise({t - (p.t_0 + p.lengthpulse + p.tend_GLY)>=0.0, 1.0 }, 0.0)
            
            deltaVt_GLY = 1.0 +  p.stim .* (p.delta_GLY .* p.KO3 .* (1.0 ./ (1.0 + exp( 1.0 .*  - p.sr_GLY .* (t - (p.t_0 + p.tstart_GLY))))) .* (1.0 - unitstepSB2)) .* gly_switch;
            Vg_glyp =  p.Vmax_glyp .* (GLYg ./ (GLYg + p.Km_GLY)) .* deltaVt_GLY;
               
            Vne_LAC =  p.Vmax_ne_LAC .* (LACn ./ (LACn + p.Km_ne_LAC) - LACe ./ (LACe + p.Km_ne_LAC));
            %Vne_LAC = max(0, Vne_LAC);

            Vge_LAC =  p.Vmax_ge_LAC .* (LACg ./ (LACg + p.Km_ge_LAC) - LACe ./ (LACe + p.Km_ge_LAC));
            Vec_LAC =  p.Vm_ec_LAC .* (LACe ./ (LACe + p.Km_ec_LAC) - LACc ./ (LACc + p.Km_ec_LAC));
            

            Veg_GLU =  p.Vmax_eg_GLU .* (GLUe ./ (GLUe + p.Km_GLU));
            Vg_gs =  p.Vmax_g_gs .* ( (GLUg ./ (GLUg + p.Km_GLU)) .* (ATPg ./ (ATPg + p.Km_ATP)));
            Vg_leak_Na =  1.06 .* (p.Sm_g ./ p.Vg) .* (p.gg_NA ./ p.F) .* ( (p.RT ./ p.F) .* log(Na_e ./ Nag) - -70.0);

            
            Vcn_O2 =  (p.PScapn ./ p.Vn) .* ( p.Ko2 .* power(p.HbOP ./ O2c - 1.0, -1.0 ./ p.nh_O2) - O2n);
            Vcg_O2 =  (p.PScapg ./ p.Vg) .* ( p.Ko2 .* power(p.HbOP ./ O2c - 1.0, -1.0 ./ p.nh_O2) - O2g);            
                                 
            Vgc_LAC =  p.Vmax_gc_LAC .* (LACg ./ (LACg + p.Km_gc_LAC) - LACc ./ (LACc + p.Km_gc_LAC));

            NADn = p.NADH_n_tot - NADHn;
            Vn_ldh =  p.kfn_ldh .* PYRn .* NADHn -  p.krn_ldh .* LACn .* NADn;

            
            NADg = p.NADH_g_tot - NADHg;
            Vg_ldh =  p.kfg_ldh .* PYRg .* NADHg -  p.krg_ldh .* LACg .* NADg;
            ADPn =  (ATPn ./ 2.0) .* ( - p.qak + power((power(p.qak, 2.0) +  4.0 .* p.qak .* (p.ATPtot ./ ATPn - 1.0)), 1.0  ./  2));
            Vn_pgk =  p.kn_pgk .* GAPn .* ADPn .* (NADn ./ NADHn);
            Vn_pk =  p.kn_pk .* PEPn .* ADPn;
            Vn_mito =  p.Vmax_n_mito .* (O2n ./ (O2n + p.Km_O2)) .* (ADPn ./ (ADPn + p.Km_ADP)) .* (PYRn ./ (PYRn + p.Km_PYR)) .* (1.0 - 1.0 ./ (1.0 + exp(  - p.aATP_mito .* ( 1.0 .* (ATPn ./ ADPn -  1.0 .* p.rATP_mito)))));
            CRn = p.PCrn_tot - PCrn;
            Vn_ck =  p.kfn_ck .* PCrn .* ADPn -  p.krn_ck .* CRn .* ATPn;
            Vn_ATPase =  p.Vmax_n_ATPase .* (ATPn ./ (ATPn + 0.0010));
            u_n = power(p.qak, 2.0) +  4.0 .* p.qak .* (p.ATPtot ./ ATPn - 1.0);
            dAMP_dATPn = (p.qak ./ 2.0 +  p.qak .* (p.ATPtot ./ ( ATPn .* power(u_n, 1.0  ./  2)))) - (1.0 +  0.50 .* power(u_n, 1.0  ./  2));
            ADPg =  (ATPg ./ 2.0) .* ( - p.qak + power((power(p.qak, 2.0) +  4.0 .* p.qak .* (p.ATPtot ./ ATPg - 1.0)), 1.0  ./  2));
            Vg_pgk =  p.kg_pgk .* GAPg .* ADPg .* (NADg./NADHg);
            Vg_pk =  p.kg_pk .* PEPg .* ADPg;
            Vg_mito =  p.Vmax_g_mito .* (O2g ./ (O2g + p.Km_O2)).*(ADPg ./ (ADPg + p.Km_ADP)).*(PYRg ./ (PYRg + p.Km_PYR)).*(1.0 - 1.0 ./ (1.0 + exp( 1.0.* - p.aATP_mito.*(ATPg ./ ADPg -  1.0.*p.rATP_mito))));

            CRg = p.PCrg_tot - PCrg;
            Vg_ck =  p.kfg_ck .* PCrg .* ADPg -  p.krg_ck .* CRg .* ATPg;
            
            Vg_ATPase =  p.Vmax_g_ATPase .* (ATPg ./ (ATPg + 0.0010));
            u_g = power(p.qak, 2.0) +  4.0 .* p.qak .* (p.ATPtot ./ ATPg - 1.0);
            dAMP_dATPg = (p.qak ./ 2.0 +  p.qak .* (p.ATPtot ./ ( ATPg .* power(u_g, 1.0  ./  2)))) - (1.0 +  0.50 .* power(u_g, 1.0  ./  2));
            
            
%             Vc_O2 =  4;
%             Vc_LAC =  ((p.R_c_cbf * CBF) ./ p.Vc) .* (p.LACa - LACc);
%             Vc_GLC =  4.64;

            % normal
            
            Vc_O2 =  ((p.R_c_cbf * CBF2) ./ p.Vc) .* (p.O2a - O2c);
            Vc_LAC =  ((p.R_c_cbf * CBF2) ./ p.Vc) .* (p.LACa - LACc);
            
            glcA = (1 - max(0, min(1.0, (t - 500)./200.0)) .* 0.5) .* p.GLCa;
            Vc_GLC =  ((p.R_c_cbf * CBF2) ./ p.Vc) .* ( p.GLCa - GLCc);
            
            % Neuron
            
            du(idx.GLCn, :) = V_en_GLC - Vn_hk;
            du(idx.G6Pn, :) = Vn_hk - Vn_pgi;
            du(idx.F6Pn, :) = Vn_pgi - Vn_pfk;
            du(idx.LACn, :) = Vn_ldh - Vne_LAC;
            du(idx.GAPn, :) = 2.0 .* Vn_pfk - Vn_pgk;
            du(idx.PYRn, :) = Vn_pk - (Vn_ldh + Vn_mito);
            du(idx.NADHn, :) = Vn_pgk - (Vn_ldh + Vn_mito);
            du(idx.O2n, :) = Vcn_O2 -  p.NAero .* Vn_mito;
            du(idx.PEPn, :) = Vn_pgk - Vn_pk;
            

            
            %Vn_pump = 0.1583;
            
            du(idx.ATPn, :) =  ((Vn_pgk + Vn_pk +  p.nOP .* Vn_mito + Vn_ck) - (Vn_hk + Vn_pfk + Vn_ATPase + Vn_pump)).* power(1.0 - dAMP_dATPn, -1.0);
            du(idx.PCrn, :) =  - Vn_ck;
            % ECS
            du(idx.GLCe, :) = Vce_GLC - (Veg_GLC .* (1.0 ./ p.Reg) +  V_en_GLC .* (1.0 ./ p.Ren));
            du(idx.LACe, :) = ( Vne_LAC .* (1.0 ./ p.Ren) +  Vge_LAC .* (1.0 ./ p.Reg)) - Vec_LAC;      
            % Astrocyte 
            
            du(idx.Nag, :) = (Vg_leak_Na +  3.0 .* Veg_GLU) -  3.0 .* Vk_pump;
            du(idx.GLCg, :) = (Vcg_GLC + Veg_GLC) - Vg_hk;
            du(idx.F6Pg, :) = Vg_pgi - Vg_pfk;
            du(idx.G6Pg, :) = (Vg_hk + Vg_glyp) - (Vg_pgi + Vg_glys);
            du(idx.GLYg, :) = Vg_glys - Vg_glyp;
            du(idx.PEPg, :) = Vg_pgk - Vg_pk;
            
            du(idx.PYRg, :) = Vg_pk - (Vg_ldh + Vg_mito);
            du(idx.NADHg, :) = Vg_pgk - (Vg_ldh + Vg_mito);
            du(idx.O2g, :) = Vcg_O2 -  p.NAero .* Vg_mito;
            du(idx.PCrg, :) =  - Vg_ck;
            du(idx.LACg, :) = Vg_ldh - (Vge_LAC + Vgc_LAC);
            du(idx.GAPg, :) =  2.0 .* Vg_pfk - Vg_pgk;
            
            du(idx.GLUg, :) = Veg_GLU - Vg_gs;
            
            du(idx.ATPg, :) = ((Vg_pgk + Vg_pk +  p.nOP .* Vg_mito + Vg_ck) - (Vg_hk + Vg_pfk + Vg_ATPase + Vk_pump + Vg_gs)) .* power(1.0 - dAMP_dATPg, -1.0);              %todo add na/K pump for astro...
            % Capillary
            du(idx.O2c, :) = Vc_O2 - ( Vcn_O2 .* (1.0 ./ p.Rcn) +  Vcg_O2 .* (1.0 ./ p.Rcg));
            du(idx.GLCc, :) = Vc_GLC - ( Vce_GLC .* (1.0 / p.Rce) +  Vcg_GLC .* (1.0 / p.Rcg));
            du(idx.LACc, :) = Vc_LAC + ( Vec_LAC .* (1.0 ./ p.Rce) +  Vgc_LAC .* (1.0 ./ p.Rcg));
            
            
            
            du = bsxfun(@times, self.enabled, du);
            
            if nargout == 2
                Uout = zeros(self.n_out, size(u, 2));

                Uout(self.idx_out.Vn_pump, :) = Vn_pump;  
             
            
            
                Uout(self.idx_out.Vk_pump, :) = Vk_pump;
                Uout(self.idx_out.deltaVt_GLY, :) = deltaVt_GLY;
                
                Uout(self.idx_out.V_en_GLC, :) = V_en_GLC;
                Uout(self.idx_out.Vn_hk, :) = Vn_hk;
                
                
                Uout(self.idx_out.Vc_O2, :) = Vc_O2;
                Uout(self.idx_out.Vc_LAC, :) = Vc_LAC;
                Uout(self.idx_out.Vc_GLC, :) = Vc_GLC;
                Uout(self.idx_out.I_pump, :) = I_pump;
                
                Uout(self.idx_out.dAMP_dATPn, :) = dAMP_dATPn;
                Uout(self.idx_out.Vn_ATPase, :) = Vn_ATPase;
                Uout(self.idx_out.Vn_pfk, :) = Vn_pfk;
                Uout(self.idx_out.Vn_hk, :) = Vn_hk;
                Uout(self.idx_out.Vn_ck, :) = Vn_ck;
                Uout(self.idx_out.Vn_mito, :) = Vn_mito;
                Uout(self.idx_out.Vn_pgk, :) = Vn_pgk;
                Uout(self.idx_out.Vn_pk, :) = Vn_pk;   
                Uout(self.idx_out.ADPn, :) = ADPn;   
                Uout(self.idx_out.Vn_ldh, :) = Vn_ldh; 
                
                Uout(self.idx_out.Vge_LAC, :) = Vge_LAC; 
                Uout(self.idx_out.Ik_pump, :) = Ik_pump; 
                Uout(self.idx_out.Ik_pump, :) = Ik_pump; 
                
                Uout(self.idx_out.CBF2, :) = CBF2; 
                
                Uout(self.idx_out.Vne_LAC, :) = Vne_LAC;
                Uout(self.idx_out.Vec_LAC, :) = Vec_LAC;
                Uout(self.idx_out.Vgc_LAC, :) = Vgc_LAC;
                Uout(self.idx_out.atp_term, :) = atp_term;
                Uout(self.idx_out.I_pumpv2, :) = I_pumpv2;
                Uout(self.idx_out.atp_term2, :) = atp_term2;
                

                Uout(self.idx_out.Vg_pgk, :) = Vg_pgk;                
                Uout(self.idx_out.Vg_pk, :) = Vg_pk;
                Uout(self.idx_out.Vg_mito, :) = Vg_mito;
                Uout(self.idx_out.Vg_ck, :) = Vg_ck;
                Uout(self.idx_out.Vg_hk, :) = Vg_hk;
                Uout(self.idx_out.Vg_pfk, :) = Vg_pfk;
                Uout(self.idx_out.Vg_ATPase, :) = Vg_ATPase;
                Uout(self.idx_out.Vg_gs, :) = Vg_gs;
                Uout(self.idx_out.GLUe, :) = GLUe;
                Uout(self.idx_out.Veg_GLU, :) = Veg_GLU;
                Uout(self.idx_out.Vg_leak_Na, :) = Vg_leak_Na;
                
                
                
                
               varargout{1} = Uout;
            end
        end 

        function [ATPn, ATPg, GLUg, Nag] = shared(self, ~,u)
           
            ATPn = u(self.index.ATPn, :);
            ATPg = u(self.index.ATPg, :);
            GLUg = u(self.index.GLUg, :);
            Nag = u(self.index.Nag, :);

        end  
        

        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
    end 
end

function idx = indices()
    % Index of parameters needing inital conditions 
    idx.GLCn = 1;            
    idx.G6Pn = 2;
    idx.F6Pn = 3;
    idx.GAPn = 4;
    idx.PEPn = 5;
    idx.PYRn = 6;
    idx.LACn = 7;
    idx.NADHn = 8;
    idx.ATPn = 9;
    idx.GLCg = 10;
    idx.F6Pg = 11;
    idx.G6Pg = 12;
    idx.GLYg = 13;
    idx.GLCe = 14;
    idx.LACe = 15;
    idx.O2c = 16;
    idx.GLCc = 17;
    idx.LACc = 18;
    idx.LACg = 19;
    idx.O2n = 20;
    idx.GAPg = 21;
    idx.PEPg = 22;
    idx.PYRg = 23;
    idx.NADHg = 24;
    idx.O2g = 25;
    idx.PCrg = 26;
    idx.ATPg = 27;
    idx.PCrn = 28;
    idx.GLUg = 29;
    idx.Nag = 30;
    
    
end

function [idx, n] = output_indices()
    % Index of all other output parameters
    idx.Vn_pump = 1;

                
    idx.Vk_pump = 2;
    idx.deltaVt_GLY = 3;
    
    idx.V_en_GLC = 4; 
    idx.Vn_hk = 5;
    idx.Vc_O2 = 6;
    idx.Vc_LAC = 7;
    idx.Vc_GLC = 8;
    idx.I_pump = 9;
    
    idx.dAMP_dATPn = 10;
    idx.Vn_ATPase = 11;
    idx.Vn_pfk = 12;
    idx.Vn_hk = 13;
    idx.Vn_ck = 14;
    idx.Vn_mito = 15;
    idx.Vn_pgk = 16;
    idx.Vn_pk = 17; 
    idx.ADPn = 18;
    
    idx.Vn_ldh = 19;
    
    idx.Vge_LAC = 20;
    idx.Ik_pump = 21;
    idx.Ik_pump = 22;
    idx.CBF2 = 23;
    idx.Vne_LAC = 24;
    idx.Vec_LAC = 25;
    idx.Vgc_LAC = 26;
    
    idx.atp_term = 27;
    idx.I_pumpv2 = 28;

    idx.atp_term2 = 29;
    
    idx.Vg_pgk = 30;
    idx.Vg_pgk = 31;
    
                  
    idx.Vg_pk = 32;
    idx.Vg_mito = 33;
    idx.Vg_ck = 34;
    idx.Vg_hk = 35;
    idx.Vg_pfk = 36;
    idx.Vg_ATPase = 37;
    idx.Vg_gs = 38;
    idx.GLUe = 39;
    idx.Veg_GLU = 40;
    idx.Vg_leak_Na = 41;
    n = numel(fieldnames(idx));
end

function params = parse_inputs(varargin)
    parser = inputParser();
    
    % Parameter for changing the wall mechanics rate constants, default 1.7
    

    parser.addParameter('nOP', 15.0);
    parser.addParameter('PCrn', 4.2529);
    parser.addParameter('NAero', 3.0);
    parser.addParameter('Rng', 1.8);
    parser.addParameter('Reg', 0.8);
    parser.addParameter('Ren', 0.4444444444444444);  
    parser.addParameter('Rcn', 0.01222);
    parser.addParameter('Rcg', 0.022);
    parser.addParameter('Rce', 0.0275);
    parser.addParameter('dHb', 0.0218);
    parser.addParameter('gn_NA', 0.0039);
    parser.addParameter('Sm_n', 40500);
    parser.addParameter('Vn', 0.45);
    parser.addParameter('Km_en_GLC', 5.32);
    parser.addParameter('Vm_en_GLC', 0.50417);
    parser.addParameter('Vmax_n_hk', 0.0513);
    parser.addParameter('Km_GLC', 0.105);
    parser.addParameter('G6P_inh_hk', 0.6);
    parser.addParameter('aG6P_inh_hk', 20.0);
    parser.addParameter('Vmaxf_n_pgi', 0.5);
    parser.addParameter('Vmaxr_n_pgi', 0.45);
    parser.addParameter('Km_G6P', 0.5);
    parser.addParameter('Km_F6P_pgi', 0.06);
    parser.addParameter('kn_pfk', 0.55783);
    parser.addParameter('Km_F6P_pfk', 0.18);
    parser.addParameter('Ki_ATP', 0.7595);
    parser.addParameter('nH', 4.0);
    parser.addParameter('kn_pgk', 0.4287);
    parser.addParameter('kn_pk', 28.6);
    parser.addParameter('kfn_ldh', 5.3);  %5.30
    parser.addParameter('krn_ldh', 0.1046); %0.1046
    parser.addParameter('Vmax_n_mito', 0.05557);
    parser.addParameter('Km_O2', 0.0029658);%0.0029658
    parser.addParameter('Km_ADP', 0.00107);
    parser.addParameter('Km_PYR', 0.0632); %0.0632
    parser.addParameter('rATP_mito', 20.0);
    parser.addParameter('aATP_mito', 5.0);
    parser.addParameter('Vmax_ne_LAC', 0.1978); %0.1978
    parser.addParameter('Km_ne_LAC', 0.09314);
    parser.addParameter('Vmax_n_ATPase', 0.04889);
    parser.addParameter('krn_ck', 0.015);
    parser.addParameter('kfn_ck', 0.0524681);
    parser.addParameter('nh_O2', 2.7);
    parser.addParameter('PScapn', 0.2202);
    parser.addParameter('Ko2', 0.089733);
    parser.addParameter('HbOP', 8.6);
    parser.addParameter('gg_NA', 0.00325);
    parser.addParameter('Sm_g', 10500);
    parser.addParameter('Vg', 0.25);
    parser.addParameter('Km_eg_GLC', 3.53);
    parser.addParameter('KO1', 1.0);
    parser.addParameter('Km_cg_GLC', 9.92);
    
    parser.addParameter('Vmax_g_hk', 0.050461);
    parser.addParameter('Vmaxf_g_pgi', 0.5);
    parser.addParameter('Vmaxr_g_pgi', 0.45);
    parser.addParameter('kg_pfk', 0.403);
    parser.addParameter('kg_pgk', 0.2514);
    parser.addParameter('kg_pk', 2.73);
    parser.addParameter('kfg_ldh', 6.2613);
    parser.addParameter('krg_ldh', 0.54682);
    parser.addParameter('Vmax_g_mito', 0.008454);
    parser.addParameter('Vmax_ge_LAC', 0.086124);
    parser.addParameter('Km_ge_LAC', 0.22163);
    parser.addParameter('Vmax_gc_LAC', 0.00021856);
    parser.addParameter('Km_gc_LAC', 0.12862);
    parser.addParameter('Vmax_g_ATPase', 0.035657);
    parser.addParameter('krg_ck', 0.02073);
    parser.addParameter('kfg_ck', 0.0243);
    parser.addParameter('PScapg', 0.2457);

    parser.addParameter('Vc', 0.0055);
 
    parser.addParameter('Km_ce_GLC', 8.4568);
    
%     parser.addParameter('Vm_cg_GLC', 0.00);
%     parser.addParameter('Vm_eg_GLC', 0.00);
%     parser.addParameter('Vm_ce_GLC', 0.0);
    parser.addParameter('Vm_cg_GLC', 0.0098394);
    parser.addParameter('Vm_eg_GLC', 0.038089);
    parser.addParameter('Vm_ce_GLC', 0.0489);
        

    parser.addParameter('GLCa', 4.8); %4.8
    parser.addParameter('LACa', 0.313);
    parser.addParameter('CO2a', 1.2);
    parser.addParameter('O2a', 8.34); %8.34
    
    parser.addParameter('Km_ec_LAC', 0.764818);
    parser.addParameter('Vm_ec_LAC', 0.0325);
    parser.addParameter('R_GLU_NA', 0.075);
    parser.addParameter('Km_GLU', 0.05);
    parser.addParameter('Vmax_g_gs', 0.3);
    parser.addParameter('Km_ATP', 0.01532);
    parser.addParameter('Vmax_eg_GLU', 3 * 0.0208); %0.0208
    
    parser.addParameter('Vmax_glys', 0.0001528);
    parser.addParameter('Km_G6P_glys', 0.5);
    parser.addParameter('GLY_inh', 4.2);
    parser.addParameter('aGLY_inh', 20.0);
    parser.addParameter('Vmax_glyp', 4.922e-5);
    parser.addParameter('Km_GLY', 1.0);
    parser.addParameter('stim', 1);
    parser.addParameter('sr_GLY', 4);
    parser.addParameter('t1', 2);
    parser.addParameter('delta_GLY', 62);
    parser.addParameter('KO3', 1);
    parser.addParameter('sr', 4.59186);
    parser.addParameter('deltaf', 0.42);
    parser.addParameter('Vv0', 0.0236);
    parser.addParameter('tv', 35.0);
    parser.addParameter('NADH_n_tot', 0.22);
    parser.addParameter('NADH_g_tot', 0.22);
    parser.addParameter('PCrn_tot', 5.0);
    parser.addParameter('PCrg_tot', 5.0);
    parser.addParameter('ATPtot', 2.379);
    parser.addParameter('qak', 0.92);
    parser.addParameter('k1', 2.22);
    parser.addParameter('k2', 0.46);
    parser.addParameter('k3', 0.43);
    parser.addParameter('dHb0', 0.064);
    parser.addParameter('t_n_stim', 2);
    parser.addParameter('v1_n', 0.041);
    parser.addParameter('v2_n', 2.55);
    

    parser.addParameter('tstart_GLY', 25); % s after input stim starts
    parser.addParameter('tend_GLY', 223); % s after input ends  
    
    
    parser.addParameter('R_init', 20); 
    parser.addParameter('CBF_init', 3.2e-2); 
    parser.addParameter('R_c_cbf', 0.43); %0.43
    parser.addParameter('ATP_init_n', 2.2595); % dimless 
    parser.addParameter('F', 9.65e4); %C mol^-1; Faraday's constant
    parser.addParameter('Imax', 0.013*6);  % mA/cm^2
    
    % To be overwritten from run script.
    parser.addParameter('startpulse', 100);    
    parser.addParameter('lengthpulse', 1);
    
    
    parser.addParameter('As', 1.586e-5);        % Surface area of soma [cm2]
    parser.addParameter('Vs', 2.16e-9);         % Volume of soma [cm3]
    
    parser.addParameter('K_init_e', 2.9);
    parser.addParameter('K_init_s', 2.9);
    parser.addParameter('Na_init_k', 13.06); %mM
    parser.addParameter('ATP_init_k', 2.24); 
    parser.addParameter('Imax_k', 0.013*6 * 0.08); % no idea about this one TODO (last part to get similar baseline as cloutier)
    parser.addParameter('R_k', 6e-8); % m
    
    
    parser.addParameter('K_Na_k', 10000); % uM
    
    parser.addParameter('ATPheight', 0.6); 
    parser.addParameter('ATPslope', 1.5); 
    parser.addParameter('ATPshift', 2.0);
    
    parser.addParameter('RT', 2577340);
    
    parser.parse(varargin{:});
    params = parser.Results;
    params.t_0 = params.startpulse;

end
function u0 = initial_conditions(idx)
    u0 = zeros(length(fieldnames(idx)), 1);
    % Inital estimations of parameters from experimental data
    
    % normal IC's
    
    u0(idx.GLCn) = 0.1529;
    u0(idx.G6Pn) = 0.7179;
    u0(idx.F6Pn) = 0.1075;
    u0(idx.GAPn) = 0.05624;
    u0(idx.PEPn) = 0.00367;
    u0(idx.PYRn) = 0.04823;
    
    u0(idx.LACn) = 0.5928; 
    u0(idx.NADHn) = 0.04109;
    u0(idx.ATPn) = 2.26;
    u0(idx.GLCg) = 0.006644;
    u0(idx.F6Pg) = 0.054;
    u0(idx.G6Pg) = 0.3911;
    
    u0(idx.GLYg) = 2.5;
    u0(idx.GLCe) = 0.22;
    u0(idx.LACe) = 0.6081;
    u0(idx.O2c) = 7.44;
    u0(idx.GLCc) = 4.63;
    u0(idx.LACc) = 0.3487;
    
    u0(idx.LACg) = 0.8066;
    u0(idx.O2n) = 0.1022;
    u0(idx.GAPg) = 0.0371;
    u0(idx.PEPg) = 0.00837;
    u0(idx.PYRg) = 0.1815;
    u0(idx.NADHg) = 0.06428;
    
    u0(idx.O2g) = 0.1599;
    u0(idx.PCrg) = 3.821;
    u0(idx.ATPg) = 1.794;
    u0(idx.PCrn) = 4.255;
    u0(idx.GLUg) = 0.0;
    u0(idx.Nag) = 13.37;
    
    
    


end
% 
% % Compute result of a piecewise function
% function x = piecewise(cases, default)
%     set = [0];
%     for i = 1:2:length(cases)
%         if (length(cases{i + 1}) == 1)
%             x(cases{i} & ~set,:) = cases{i + 1};
%         else
%             x(cases{i} & ~set,:) = cases{i + 1}(cases{i} & ~set);
%         end
%         set = set | cases{i};
%         if(set), break, end
%     end
%     if (length(default) == 1)
%         x(~set,:) = default;
%     else
%         x(~set,:) = default(~set);
%     end
% end
