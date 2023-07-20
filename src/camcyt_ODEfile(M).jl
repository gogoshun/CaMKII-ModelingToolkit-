# This function describes the ODE's for CaM, CaMKII, and CaN.
using DifferentialEquations
using ModelingToolkit

export get_camcyt_equations

const Mg = 1          # [mM]
const expression = 1  # expression = 1 ("WT"), =2 ("OE"), or =3 ("KO")

function get_camcyt_equations()
    @variables t
    D = Differential(t)

    @variables CaM_cyt(t) Ca2CaM_cyt(t) Ca4CaM_cyt(t) CaMB_cyt(t) Ca2CaMB_cyt(t) Ca4CaMB_cyt(t) Pb2_cyt(t) Pb_cyt(t) Pt_cyt(t) Pt2_cyt(t) 
    @variables Pa_cyt(t) Ca4CaN_cyt(t) CaMCa4CaN_cyt(t) Ca2CaMCa4CaN_cyt(t) Ca4CaMCa4CaN_cyt(t) JCaCyt(t) Ca_Cyt(t) Cai(t) 

    ## Parameters    
    @parameters K = 135         # [mM]
    @parameters Btot = 24.2     # [uM]
    @parameters CaMKIItot = 120*8.293e-4  # [uM]
    @parameters CaNtot = 3e-3             # [uM]
    @parameters PP1tot = 0.57             # [uM]
    ## ADJUST CaMKII ACTIVITY LEVELS (expression = "WT", "OE", or "KO")
    if expression == 2          # "OE"
        @parameters n_OE=6
        @parameters CaMKIItot = 120*8.293e-4*n_OE     # [uM]
    elseif expression == 3      # "KO"
        @parameters CaMKIItot = 0                     # [uM]
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
    @variables rcn02_cyt(t) rcn24_cyt(t) B_cyt(t) rcn02B_cyt(t) rcn24B_cyt(t) rcn0B_cyt(t) rcn2B_cyt(t) rcn4B_cyt(t) Ca2CaN_cyt(t) rcnCa4CaN_cyt(t) rcn02CaN_cyt(t) rcn24CaN_cyt(t) rcn0CaN_cyt(t) rcn2CaN_cyt(t) rcn4CaN_cyt(t)
    @variables Pi_cyt(t) rcnCKib2_cyt(t) rcnCKb2b_cyt(t) rcnCKib_cyt(t) T_cyt(t) kbt_cyt(t) rcnCKbt_cyt(t) rcnCKtt2_cyt(t) rcnCKta_cyt(t) rcnCKt2a_cyt(t) rcnCKt2b2_cyt(t) rcnCKai_cyt(t)
    
    fluxes_eqs = [
        Ca_Cyt ~ Cai * 1e3,
        # CaM Reaction fluxes
        rcn02_cyt ~ k02*Ca_Cyt^2*CaM_cyt - k20*Ca2CaM_cyt,
        rcn24_cyt ~ k24*Ca_Cyt^2*Ca2CaM_cyt - k42*Ca4CaM_cyt,
        # CaM buffer fluxes
        B_cyt ~ Btot - CaMB_cyt - Ca2CaMB_cyt - Ca4CaMB_cyt,
        rcn02B_cyt ~ k02B*Ca_Cyt^2*CaMB_cyt - k20B*Ca2CaMB_cyt,
        rcn24B_cyt ~ k24B*Ca_Cyt^2*Ca2CaMB_cyt - k42B*Ca4CaMB_cyt,
        rcn0B_cyt ~ k0Bon*CaM_cyt*B_cyt - k0Boff*CaMB_cyt,
        rcn2B_cyt ~ k2Bon*Ca2CaM_cyt*B_cyt - k2Boff*Ca2CaMB_cyt,
        rcn4B_cyt ~ k4Bon*Ca4CaM_cyt*B_cyt - k4Boff*Ca4CaMB_cyt,
        # CaN reaction fluxes 
        Ca2CaN_cyt ~ CaNtot - Ca4CaN_cyt - CaMCa4CaN_cyt - Ca2CaMCa4CaN_cyt - Ca4CaMCa4CaN_cyt,
        rcnCa4CaN_cyt ~ kcanCaon*Ca_Cyt^2*Ca2CaN_cyt - kcanCaoff*Ca4CaN_cyt,
        rcn02CaN_cyt ~ k02can*Ca_Cyt^2*CaMCa4CaN_cyt - k20can*Ca2CaMCa4CaN_cyt,
        rcn24CaN_cyt ~ k24can*Ca_Cyt^2*Ca2CaMCa4CaN_cyt - k42can*Ca4CaMCa4CaN_cyt,
        rcn0CaN_cyt ~ kcanCaM0on*CaM_cyt*Ca4CaN_cyt - kcanCaM0off*CaMCa4CaN_cyt,
        rcn2CaN_cyt ~ kcanCaM2on*Ca2CaM_cyt*Ca4CaN_cyt - kcanCaM2off*Ca2CaMCa4CaN_cyt,
        rcn4CaN_cyt ~ kcanCaM4on*Ca4CaM_cyt*Ca4CaN_cyt - kcanCaM4off*Ca4CaMCa4CaN_cyt,
        # CaMKII reaction fluxes
        Pi_cyt ~ 1 - Pb2_cyt - Pb_cyt - Pt_cyt - Pt2_cyt - Pa_cyt,
        rcnCKib2_cyt ~ kib2*Ca2CaM_cyt*Pi_cyt - kb2i*Pb2_cyt,
        rcnCKb2b_cyt ~ kb24*Ca_Cyt^2*Pb2_cyt - kb42*Pb_cyt,
        rcnCKib_cyt ~ kib*Ca4CaM_cyt*Pi_cyt - kbi*Pb_cyt,
        T_cyt ~ Pb_cyt + Pt_cyt + Pt2_cyt + Pa_cyt,
        kbt_cyt ~ 0.055*T_cyt + 0.0074*T_cyt^2 + 0.015*T_cyt^3,
        rcnCKbt_cyt ~ kbt_cyt*Pb_cyt - kpp1*PP1tot*Pt_cyt/(Kmpp1+CaMKIItot*Pt_cyt),
        rcnCKtt2_cyt ~ kt42*Pt_cyt - kt24*Ca_Cyt^2*Pt2_cyt,
        rcnCKta_cyt ~ kta*Pt_cyt - kat*Ca4CaM_cyt*Pa_cyt,
        rcnCKt2a_cyt ~ kt2a*Pt2_cyt - kat2*Ca2CaM_cyt*Pa_cyt,
        rcnCKt2b2_cyt ~ kpp1*PP1tot*Pt2_cyt/(Kmpp1+CaMKIItot*Pt2_cyt),
        rcnCKai_cyt ~ kpp1*PP1tot*Pa_cyt/(Kmpp1+CaMKIItot*Pa_cyt)
    ]
    

    ## Ordinary Differential Equations
    @parameters Vmyo = 2.1454e-11           # [L]
    @parameters kSLmyo = 8.587e-15          # [L/msec]

    @variables CaM_sl(t) Ca2CaM_sl(t) Ca4CaM_sl(t) J_cam_SLmyo(t) J_ca2cam_SLmyo(t) J_ca4cam_SLmyo(t)

    # CaM equations
    CaMSL_eqs = [
        J_cam_SLmyo ~ kSLmyo*(CaM_sl-CaM_cyt),                                          # [umol/msec]
        J_ca2cam_SLmyo ~ kSLmyo*(Ca2CaM_sl-Ca2CaM_cyt),                                 # [umol/msec]
        J_ca4cam_SLmyo ~ kSLmyo*(Ca4CaM_sl-Ca4CaM_cyt),                                 # [umol/msec]
        D(CaM_cyt) ~ 1e-3*(-rcn02_cyt - rcn0B_cyt - rcn0CaN_cyt) + J_cam_SLmyo/Vmyo,                                                                    # du[1]
        D(Ca2CaM_cyt) ~ 1e-3*(rcn02_cyt - rcn24_cyt - rcn2B_cyt - rcn2CaN_cyt + CaMKIItot.*(-rcnCKib2_cyt + rcnCKt2a_cyt)) + J_ca2cam_SLmyo/Vmyo,       # du[2]
        D(Ca4CaM_cyt) ~ 1e-3*(rcn24_cyt - rcn4B_cyt - rcn4CaN_cyt + CaMKIItot.*(-rcnCKib_cyt+rcnCKta_cyt)) + J_ca4cam_SLmyo/Vmyo,                       # du[3]
        D(CaMB_cyt) ~ 1e-3*(rcn0B_cyt-rcn02B_cyt),                      # du[4]
        D(Ca2CaMB_cyt) ~ 1e-3*(rcn02B_cyt + rcn2B_cyt - rcn24B_cyt),    # du[5]
        D(Ca4CaMB_cyt) ~ 1e-3*(rcn24B_cyt + rcn4B_cyt),                 # du[6]
        # CaMKII equations
        D(Pb2_cyt) ~ 1e-3*(rcnCKib2_cyt - rcnCKb2b_cyt + rcnCKt2b2_cyt),     # du[7]
        D(Pb_cyt) ~ 1e-3*(rcnCKib_cyt + rcnCKb2b_cyt - rcnCKbt_cyt),         # du[8]
        D(Pt_cyt) ~ 1e-3*(rcnCKbt_cyt-rcnCKta_cyt-rcnCKtt2_cyt),             # du[9]
        D(Pt2_cyt) ~ 1e-3*(rcnCKtt2_cyt-rcnCKt2a_cyt-rcnCKt2b2_cyt),         # du[10]
        D(Pa_cyt) ~ 1e-3*(rcnCKta_cyt+rcnCKt2a_cyt-rcnCKai_cyt),             # du[11]
        # CaN equations
        D(Ca4CaN_cyt) ~ 1e-3*(rcnCa4CaN_cyt - rcn0CaN_cyt - rcn2CaN_cyt - rcn4CaN_cyt),      # du[12]
        D(CaMCa4CaN_cyt) ~ 1e-3*(rcn0CaN_cyt - rcn02CaN_cyt),                        # du[13]
        D(Ca2CaMCa4CaN_cyt) ~ 1e-3*(rcn2CaN_cyt+rcn02CaN_cyt-rcn24CaN_cyt),              # du[14]
        D(Ca4CaMCa4CaN_cyt) ~ 1e-3*(rcn4CaN_cyt+rcn24CaN_cyt),                       # du[15]
        ## For adjusting Ca buffering in EC coupling model
        JCaCyt ~ 1e-3*(2*CaMKIItot*(rcnCKtt2_cyt-rcnCKb2b_cyt) - 2*(rcn02_cyt+rcn24_cyt+rcn02B_cyt+rcn24B_cyt+rcnCa4CaN_cyt+rcn02CaN_cyt+rcn24CaN_cyt))   # [uM/msec]
    ]
    
    eqs = append!(fluxes_eqs, CaMSL_eqs)

    @named sys = ODESystem(eqs, t)

    return sys
end
