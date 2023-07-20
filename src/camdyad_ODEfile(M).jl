# This function describes the ODE's for CaM, CaMKII, and CaN.
using DifferentialEquations
using ModelingToolkit

export get_camdyad_equations

const Mg = 1      # [mM]
const expression = 1  # expression = 1 ("WT"), =2 ("OE"), or =3 ("KO")
function get_camdyad_equations()
    @variables t
    D = Differential(t)

    @variables CaM_dyad(t) Ca2CaM_dyad(t) Ca4CaM_dyad(t) CaMB_dyad(t) Ca2CaMB_dyad(t) Ca4CaMB_dyad(t) Pb2_dyad(t) Pb_dyad(t) Pt_dyad(t) Pt2_dyad(t) 
    @variables Pa_dyad(t) Ca4CaN_dyad(t) CaMCa4CaN_dyad(t) Ca2CaMCa4CaN_dyad(t) Ca4CaMCa4CaN_dyad(t) JCaDyad(t) Ca_Dyad(t) 
    @variables Ca_j(t)

    ## Parameters
    
    @parameters K = 135     # [mM]
    @parameters Btot = 0
    @parameters CaMKIItot = 120         # [uM]
    @parameters CaNtot = 3e-3/8.293e-4  # [uM]
    @parameters PP1tot = 96.5           # [uM]

    ## ADJUST CaMKII ACTIVITY LEVELS (expression = "WT", "OE", or "KO")
    
    
    if expression == 2          # "OE"
        @parameters n_OE=6
        @parameters CaMKIItot = 120*n_OE            # [uM] 
    elseif expression == 3      # "KO"
        @parameters CaMKIItot = 0                   # [uM] 
    end

    ## Parameters
    # Ca/CaM parameters
    if Mg <= 1
        @parameters Kd02 = 0.0025*(1+K/0.94-Mg/0.012)*(1+K/8.1+Mg/0.022)  # [uM^2]
        @parameters Kd24 = 0.128*(1+K/0.64+Mg/0.0014)*(1+K/13.0-Mg/0.153) # [uM^2]
    else
        @parameters Kd02 = 0.0025*(1+K/0.94-1/0.012+(Mg-1)/0.060)*(1+K/8.1+1/0.022+(Mg-1)/0.068)   # [uM^2]
        @parameters Kd24 = 0.128*(1+K/0.64+1/0.0014+(Mg-1)/0.005)*(1+K/13.0-1/0.153+(Mg-1)/0.150)  # [uM^2]
    end
    @parameters k20 = 10               # [s^-1]      
    @parameters k02 = k20/Kd02         # [uM^-2 s^-1]
    @parameters k42 = 500              # [s^-1]      
    @parameters k24 = k42/Kd24         # [uM^-2 s^-1]
    
    # CaM buffering (B) parameters
    @parameters k0Boff = 0.0014        # [s^-1] 
    @parameters k0Bon = k0Boff/0.2   # [uM^-1 s^-1] kon = koff/Kd
    @parameters k2Boff = k0Boff/100    # [s^-1] 
    @parameters k2Bon = k0Bon          # [uM^-1 s^-1]
    @parameters k4Boff = k2Boff        # [s^-1]
    @parameters k4Bon = k0Bon          # [uM^-1 s^-1]

    # using thermodynamic constraints
    @parameters k20B = k20/100 # [s^-1] thermo constraint on loop 1
    @parameters k02B = k02     # [uM^-2 s^-1] 
    @parameters k42B = k42     # [s^-1] thermo constraint on loop 2
    @parameters k24B = k24     # [uM^-2 s^-1]
    
    # CaMKII parameters
    # Wi Wa Wt Wp
    @parameters kbi = 2.2      # [s^-1] (Ca4CaM dissocation from Wb)
    @parameters kib = kbi/33.5e-3 # [uM^-1 s^-1]
    @parameters kib2 = kib
    @parameters kb2i = kib2*5
    @parameters kb24 = k24
    @parameters kb42 = k42*33.5e-3/5
    @parameters kpp1 = 1.72     # [s^-1] (PP1-dep dephosphorylation rates)
    @parameters Kmpp1 = 11.5    # [uM]
    @parameters kta = kbi/1000  # [s^-1] (Ca4CaM dissociation from Wt)
    @parameters kat = kib       # [uM^-1 s^-1] (Ca4CaM reassociation with Wa)
    @parameters kt42 = k42*33.5e-6/5
    @parameters kt24 = k24
    @parameters kat2 = kib
    @parameters kt2a = kib*5
    
    # CaN parameters
    @parameters kcanCaoff = 1              # [s^-1] 
    @parameters kcanCaon = kcanCaoff/0.5   # [uM^-1 s^-1] 
    @parameters kcanCaM4on = 46            # [uM^-1 s^-1]
    @parameters kcanCaM4off = 1.3e-3       # [s^-1]
    @parameters kcanCaM2on = kcanCaM4on
    @parameters kcanCaM2off = 2508*kcanCaM4off
    @parameters kcanCaM0on = kcanCaM4on
    @parameters kcanCaM0off = 165*kcanCaM2off
    @parameters k02can = k02
    @parameters k20can = k20/165
    @parameters k24can = k24
    @parameters k42can = k20/2508
 
    ## Fluxes
    @variables rcn02_dyad(t) rcn24_dyad(t) B_dyad(t) rcn02B_dyad(t) rcn24B_dyad(t) rcn0B_dyad(t) rcn2B_dyad(t) rcn4B_dyad(t) Ca2CaN_dyad(t) rcnCa4CaN_dyad(t) rcn02CaN_dyad(t) rcn24CaN_dyad(t) rcn0CaN_dyad(t) rcn2CaN_dyad(t) rcn4CaN_dyad(t)
    @variables Pi_dyad(t) rcnCKib2_dyad(t) rcnCKb2b_dyad(t) rcnCKib_dyad(t) T_dyad(t) kbt_dyad(t) rcnCKbt_dyad(t) rcnCKtt2_dyad(t) rcnCKta_dyad(t) rcnCKt2a_dyad(t) rcnCKt2b2_dyad(t) rcnCKai_dyad(t)
    
    fluxes_eqs = [
        Ca_Dyad ~ Ca_j * 1e3,
        # CaM Reaction fluxes
        rcn02_dyad ~ k02*Ca_Dyad^2*CaM_dyad - k20*Ca2CaM_dyad,
        rcn24_dyad ~ k24*Ca_Dyad^2*Ca2CaM_dyad - k42*Ca4CaM_dyad,
        # CaM buffer fluxes
        B_dyad ~ Btot - CaMB_dyad - Ca2CaMB_dyad - Ca4CaMB_dyad,
        rcn02B_dyad ~ k02B*Ca_Dyad^2*CaMB_dyad - k20B*Ca2CaMB_dyad,
        rcn24B_dyad ~ k24B*Ca_Dyad^2*Ca2CaMB_dyad - k42B*Ca4CaMB_dyad,
        rcn0B_dyad ~ k0Bon*CaM_dyad*B_dyad - k0Boff*CaMB_dyad,
        rcn2B_dyad ~ k2Bon*Ca2CaM_dyad*B_dyad - k2Boff*Ca2CaMB_dyad,
        rcn4B_dyad ~ k4Bon*Ca4CaM_dyad*B_dyad - k4Boff*Ca4CaMB_dyad,
        # CaN reaction fluxes 
        Ca2CaN_dyad ~ CaNtot - Ca4CaN_dyad - CaMCa4CaN_dyad - Ca2CaMCa4CaN_dyad - Ca4CaMCa4CaN_dyad,
        rcnCa4CaN_dyad ~ kcanCaon*Ca_Dyad^2*Ca2CaN_dyad - kcanCaoff*Ca4CaN_dyad,
        rcn02CaN_dyad ~ k02can*Ca_Dyad^2*CaMCa4CaN_dyad - k20can*Ca2CaMCa4CaN_dyad,
        rcn24CaN_dyad ~ k24can*Ca_Dyad^2*Ca2CaMCa4CaN_dyad - k42can*Ca4CaMCa4CaN_dyad,
        rcn0CaN_dyad ~ kcanCaM0on*CaM_dyad*Ca4CaN_dyad - kcanCaM0off*CaMCa4CaN_dyad,
        rcn2CaN_dyad ~ kcanCaM2on*Ca2CaM_dyad*Ca4CaN_dyad - kcanCaM2off*Ca2CaMCa4CaN_dyad,
        rcn4CaN_dyad ~ kcanCaM4on*Ca4CaM_dyad*Ca4CaN_dyad - kcanCaM4off*Ca4CaMCa4CaN_dyad,
        # CaMKII reaction fluxes
        Pi_dyad ~ 1 - Pb2_dyad - Pb_dyad - Pt_dyad - Pt2_dyad - Pa_dyad,
        rcnCKib2_dyad ~ kib2*Ca2CaM_dyad*Pi_dyad - kb2i*Pb2_dyad,
        rcnCKb2b_dyad ~ kb24*Ca_Dyad^2*Pb2_dyad - kb42*Pb_dyad,
        rcnCKib_dyad ~ kib*Ca4CaM_dyad*Pi_dyad - kbi*Pb_dyad,
        T_dyad ~ Pb_dyad + Pt_dyad + Pt2_dyad + Pa_dyad,
        kbt_dyad ~ 0.055*T_dyad + 0.0074*T_dyad^2 + 0.015*T_dyad^3,
        rcnCKbt_dyad ~ kbt_dyad*Pb_dyad - kpp1*PP1tot*Pt_dyad/(Kmpp1+CaMKIItot*Pt_dyad),
        rcnCKtt2_dyad ~ kt42*Pt_dyad - kt24*Ca_Dyad^2*Pt2_dyad,
        rcnCKta_dyad ~ kta*Pt_dyad - kat*Ca4CaM_dyad*Pa_dyad,
        rcnCKt2a_dyad ~ kt2a*Pt2_dyad - kat2*Ca2CaM_dyad*Pa_dyad,
        rcnCKt2b2_dyad ~ kpp1*PP1tot*Pt2_dyad/(Kmpp1+CaMKIItot*Pt2_dyad),
        rcnCKai_dyad ~ kpp1*PP1tot*Pa_dyad/(Kmpp1+CaMKIItot*Pa_dyad)
    ]
    
    ## Ordinary Differential Equations
    @parameters Vmyo = 2.1454e-11           # [L]
    @parameters Vdyad = 1.7790e-014         # [L]
    @parameters VSL = 6.6013e-013           # [L]
    @parameters kDyadSL = 3.6363e-16 	    # [L/msec]
    @parameters kSLmyo = 8.587e-15          # [L/msec]
    @parameters CaMKIItotDyad = 120         # [uM]
    @parameters BtotDyad = 1.54/8.293e-4    # [uM]

    @variables CaMtotDyad(t) Bdyad(t) J_cam_dyadSL(t) J_ca2cam_dyadSL(t) J_ca4cam_dyadSL(t) CaM_sl(t) Ca2CaM_sl(t) Ca4CaM_sl(t)


    # CaM equations
    CaMSL_eqs = [
        CaMtotDyad ~ CaM_dyad+Ca2CaM_dyad+Ca4CaM_dyad+CaMB_dyad+Ca2CaMB_dyad+Ca4CaMB_dyad+CaMKIItotDyad*(Pb2_dyad+Pb_dyad+Pt_dyad+Pt2_dyad)+CaMCa4CaN_dyad+Ca2CaMCa4CaN_dyad+Ca4CaMCa4CaN_dyad,
        Bdyad ~ BtotDyad - CaMtotDyad,                                                  # [uM dyad]
        J_cam_dyadSL ~ 1e-3*(k0Boff*CaM_dyad - k0Bon*Bdyad*CaM_sl),                     # [uM/msec dyad]
        J_ca2cam_dyadSL ~ 1e-3*(k2Boff*Ca2CaM_dyad - k2Bon*Bdyad*Ca2CaM_sl),            # [uM/msec dyad]
        J_ca4cam_dyadSL ~ 1e-3*(k2Boff*Ca4CaM_dyad - k4Bon*Bdyad*Ca4CaM_sl),            # [uM/msec dyad]
        D(CaM_dyad) ~ 1e-3*(-rcn02_dyad - rcn0B_dyad - rcn0CaN_dyad)-J_cam_dyadSL,                                                                      # du[1]
        D(Ca2CaM_dyad) ~ 1e-3*(rcn02_dyad - rcn24_dyad - rcn2B_dyad - rcn2CaN_dyad + CaMKIItot.*(-rcnCKib2_dyad + rcnCKt2a_dyad))-J_ca2cam_dyadSL,      # du[2]
        D(Ca4CaM_dyad) ~ 1e-3*(rcn24_dyad - rcn4B_dyad - rcn4CaN_dyad + CaMKIItot.*(-rcnCKib_dyad+rcnCKta_dyad))-J_ca4cam_dyadSL,                       # du[3]
        D(CaMB_dyad) ~ 1e-3*(rcn0B_dyad-rcn02B_dyad),                       # du[4]
        D(Ca2CaMB_dyad) ~ 1e-3*(rcn02B_dyad + rcn2B_dyad - rcn24B_dyad),    # du[5]
        D(Ca4CaMB_dyad) ~ 1e-3*(rcn24B_dyad + rcn4B_dyad),                  # du[6]
        # CaMKII equations
        D(Pb2_dyad) ~ 1e-3*(rcnCKib2_dyad - rcnCKb2b_dyad + rcnCKt2b2_dyad),     # du[7]
        D(Pb_dyad) ~ 1e-3*(rcnCKib_dyad + rcnCKb2b_dyad - rcnCKbt_dyad),         # du[8]
        D(Pt_dyad) ~ 1e-3*(rcnCKbt_dyad-rcnCKta_dyad-rcnCKtt2_dyad),             # du[9]
        D(Pt2_dyad) ~ 1e-3*(rcnCKtt2_dyad-rcnCKt2a_dyad-rcnCKt2b2_dyad),         # du[10]
        D(Pa_dyad) ~ 1e-3*(rcnCKta_dyad+rcnCKt2a_dyad-rcnCKai_dyad),             # du[11]
        # CaN equations
        D(Ca4CaN_dyad) ~ 1e-3*(rcnCa4CaN_dyad - rcn0CaN_dyad - rcn2CaN_dyad - rcn4CaN_dyad),      # du[12]
        D(CaMCa4CaN_dyad) ~ 1e-3*(rcn0CaN_dyad - rcn02CaN_dyad),                        # du[13]
        D(Ca2CaMCa4CaN_dyad) ~ 1e-3*(rcn2CaN_dyad+rcn02CaN_dyad-rcn24CaN_dyad),              # du[14]
        D(Ca4CaMCa4CaN_dyad) ~ 1e-3*(rcn4CaN_dyad+rcn24CaN_dyad),                       # du[15]
        ## For adjusting Ca buffering in EC coupling model
        JCaDyad ~ 1e-3*(2*CaMKIItot*(rcnCKtt2_dyad-rcnCKb2b_dyad) - 2*(rcn02_dyad+rcn24_dyad+rcn02B_dyad+rcn24B_dyad+rcnCa4CaN_dyad+rcn02CaN_dyad+rcn24CaN_dyad))   # [uM/msec]
    ]
    
    eqs = append!(fluxes_eqs, CaMSL_eqs)

    @named sys = ODESystem(eqs, t)

    return sys
end
