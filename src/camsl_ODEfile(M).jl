# This function describes the ODE's for CaM, CaMKII, and CaN.
using DifferentialEquations
using ModelingToolkit

export get_camsl_equations

const expression = 1        # expression = 1 ("WT"), =2 ("OE"), or =3 ("KO")
const Mg = 1                # [mM]

function get_camsl_equations()
    @variables t
    D = Differential(t)

    @variables CaM_sl(t) Ca2CaM_sl(t) Ca4CaM_sl(t) CaMB_sl(t) Ca2CaMB_sl(t) Ca4CaMB_sl(t) Pb2_sl(t) Pb_sl(t) Pt_sl(t) Pt2_sl(t) 
    @variables Pa_sl(t) Ca4CaN_sl(t) CaMCa4CaN_sl(t) Ca2CaMCa4CaN_sl(t) Ca4CaMCa4CaN_sl(t) JCaSL(t) Ca_SL(t) Ca_sl(t)

    ## Parameters    
    @parameters K = 135         # [mM]
    @parameters Btot = 24.2     # [uM]
    @parameters CaMKIItot = 120*8.293e-4  # [uM]
    @parameters CaNtot = 3e-3             # [uM]
    @parameters PP1tot = 0.57             # [uM]

    ## ADJUST CaMKII ACTIVITY LEVELS (expression = "WT", "OE", or "KO")
    if expression == 2          # "OE"
        n_OE=6
        CaMKIItot = 120*8.293e-4*n_OE     # [uM]
    elseif expression == 3      # "KO"
        CaMKIItot = 0                     # [uM]
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
    @variables rcn02_sl(t) rcn24_sl(t) B_sl(t) rcn02B_sl(t) rcn24B_sl(t) rcn0B_sl(t) rcn2B_sl(t) rcn4B_sl(t) Ca2CaN_sl(t) rcnCa4CaN_sl(t) rcn02CaN_sl(t) rcn24CaN_sl(t) rcn0CaN_sl(t) rcn2CaN_sl(t) rcn4CaN_sl(t)
    @variables Pi_sl(t) rcnCKib2_sl(t) rcnCKb2b_sl(t) rcnCKib_sl(t) T_sl(t) kbt_sl(t) rcnCKbt_sl(t) rcnCKtt2_sl(t) rcnCKta_sl(t) rcnCKt2a_sl(t) rcnCKt2b2_sl(t) rcnCKai_sl(t)
    
    fluxes_eqs = [
        Ca_SL ~ Ca_sl * 1e3,
        # CaM Reaction fluxes
        rcn02_sl ~ k02*Ca_SL^2*CaM_sl - k20*Ca2CaM_sl,
        rcn24_sl ~ k24*Ca_SL^2*Ca2CaM_sl - k42*Ca4CaM_sl,
        # CaM buffer fluxes
        B_sl ~ Btot - CaMB_sl - Ca2CaMB_sl - Ca4CaMB_sl,
        rcn02B_sl ~ k02B*Ca_SL^2*CaMB_sl - k20B*Ca2CaMB_sl,
        rcn24B_sl ~ k24B*Ca_SL^2*Ca2CaMB_sl - k42B*Ca4CaMB_sl,
        rcn0B_sl ~ k0Bon*CaM_sl*B_sl - k0Boff*CaMB_sl,
        rcn2B_sl ~ k2Bon*Ca2CaM_sl*B_sl - k2Boff*Ca2CaMB_sl,
        rcn4B_sl ~ k4Bon*Ca4CaM_sl*B_sl - k4Boff*Ca4CaMB_sl,
        # CaN reaction fluxes 
        Ca2CaN_sl ~ CaNtot - Ca4CaN_sl - CaMCa4CaN_sl - Ca2CaMCa4CaN_sl - Ca4CaMCa4CaN_sl,
        rcnCa4CaN_sl ~ kcanCaon*Ca_SL^2*Ca2CaN_sl - kcanCaoff*Ca4CaN_sl,
        rcn02CaN_sl ~ k02can*Ca_SL^2*CaMCa4CaN_sl - k20can*Ca2CaMCa4CaN_sl,
        rcn24CaN_sl ~ k24can*Ca_SL^2*Ca2CaMCa4CaN_sl - k42can*Ca4CaMCa4CaN_sl,
        rcn0CaN_sl ~ kcanCaM0on*CaM_sl*Ca4CaN_sl - kcanCaM0off*CaMCa4CaN_sl,
        rcn2CaN_sl ~ kcanCaM2on*Ca2CaM_sl*Ca4CaN_sl - kcanCaM2off*Ca2CaMCa4CaN_sl,
        rcn4CaN_sl ~ kcanCaM4on*Ca4CaM_sl*Ca4CaN_sl - kcanCaM4off*Ca4CaMCa4CaN_sl,
        # CaMKII reaction fluxes
        Pi_sl ~ 1 - Pb2_sl - Pb_sl - Pt_sl - Pt2_sl - Pa_sl,
        rcnCKib2_sl ~ kib2*Ca2CaM_sl*Pi_sl - kb2i*Pb2_sl,
        rcnCKb2b_sl ~ kb24*Ca_SL^2*Pb2_sl - kb42*Pb_sl,
        rcnCKib_sl ~ kib*Ca4CaM_sl*Pi_sl - kbi*Pb_sl,
        T_sl ~ Pb_sl + Pt_sl + Pt2_sl + Pa_sl,
        kbt_sl ~ 0.055*T_sl + 0.0074*T_sl^2 + 0.015*T_sl^3,
        rcnCKbt_sl ~ kbt_sl*Pb_sl - kpp1*PP1tot*Pt_sl/(Kmpp1+CaMKIItot*Pt_sl),
        rcnCKtt2_sl ~ kt42*Pt_sl - kt24*Ca_SL^2*Pt2_sl,
        rcnCKta_sl ~ kta*Pt_sl - kat*Ca4CaM_sl*Pa_sl,
        rcnCKt2a_sl ~ kt2a*Pt2_sl - kat2*Ca2CaM_sl*Pa_sl,
        rcnCKt2b2_sl ~ kpp1*PP1tot*Pt2_sl/(Kmpp1+CaMKIItot*Pt2_sl),
        rcnCKai_sl ~ kpp1*PP1tot*Pa_sl/(Kmpp1+CaMKIItot*Pa_sl)
    ]
    

    ## Ordinary Differential Equations
    @parameters Vmyo = 2.1454e-11           # [L]
    @parameters Vdyad = 1.7790e-014         # [L]
    @parameters VSL = 6.6013e-013           # [L]
    @parameters kDyadSL = 3.6363e-16 	    # [L/msec]
    @parameters kSLmyo = 8.587e-15          # [L/msec]
    @parameters CaMKIItotDyad = 120         # [uM]
    @parameters BtotDyad = 1.54/8.293e-4    # [uM]
 
    @variables CaMtotDyad(t) Bdyad(t) J_cam_dyadSL(t) J_ca2cam_dyadSL(t) J_ca4cam_dyadSL(t) J_cam_SLmyo(t) J_ca2cam_SLmyo(t) J_ca4cam_SLmyo(t)
    @variables CaM_dyad(t) Ca2CaM_dyad(t) Ca4CaM_dyad(t) CaMB_dyad(t) Ca2CaMB_dyad(t) Ca4CaMB_dyad(t) Pb2_dyad(t) Pb_dyad(t) Pt_dyad(t) Pt2_dyad(t) 
    @variables CaMCa4CaN_dyad(t) Ca2CaMCa4CaN_dyad(t) Ca4CaMCa4CaN_dyad(t) CaM_cyt(t) Ca2CaM_cyt(t) Ca4CaM_cyt(t)


    # CaM equations
    CaMSL_eqs = [
        CaMtotDyad ~ CaM_dyad+Ca2CaM_dyad+Ca4CaM_dyad+CaMB_dyad+Ca2CaMB_dyad+Ca4CaMB_dyad+CaMKIItotDyad*(Pb2_dyad+Pb_dyad+Pt_dyad+Pt2_dyad)+CaMCa4CaN_dyad+Ca2CaMCa4CaN_dyad+Ca4CaMCa4CaN_dyad,
        Bdyad ~ BtotDyad - CaMtotDyad,                                                  # [uM dyad]
        J_cam_dyadSL ~ 1e-3*(k0Boff*CaM_dyad - k0Bon*Bdyad*CaM_sl),                     # [uM/msec dyad]
        J_ca2cam_dyadSL ~ 1e-3*(k2Boff*Ca2CaM_dyad - k2Bon*Bdyad*Ca2CaM_sl),            # [uM/msec dyad]
        J_ca4cam_dyadSL ~ 1e-3*(k2Boff*Ca4CaM_dyad - k4Bon*Bdyad*Ca4CaM_sl),            # [uM/msec dyad]
        J_cam_SLmyo ~ kSLmyo*(CaM_sl-CaM_cyt),                                          # [umol/msec]
        J_ca2cam_SLmyo ~ kSLmyo*(Ca2CaM_sl-Ca2CaM_cyt),                                 # [umol/msec]
        J_ca4cam_SLmyo ~ kSLmyo*(Ca4CaM_sl-Ca4CaM_cyt),                                 # [umol/msec]
        D(CaM_sl) ~ 1e-3*(-rcn02_sl - rcn0B_sl - rcn0CaN_sl)+ J_cam_dyadSL*Vdyad/VSL - J_cam_SLmyo/VSL,                                                                     # du[1]
        D(Ca2CaM_sl) ~ 1e-3*(rcn02_sl - rcn24_sl - rcn2B_sl - rcn2CaN_sl + CaMKIItot.*(-rcnCKib2_sl + rcnCKt2a_sl))+ J_ca2cam_dyadSL*Vdyad/VSL - J_ca2cam_SLmyo/VSL,        # du[2]
        D(Ca4CaM_sl) ~ 1e-3*(rcn24_sl - rcn4B_sl - rcn4CaN_sl + CaMKIItot.*(-rcnCKib_sl+rcnCKta_sl))+ J_ca4cam_dyadSL*Vdyad/VSL - J_ca4cam_SLmyo/VSL,                       # du[3]
        D(CaMB_sl) ~ 1e-3*(rcn0B_sl-rcn02B_sl),                     # du[4]
        D(Ca2CaMB_sl) ~ 1e-3*(rcn02B_sl + rcn2B_sl - rcn24B_sl),    # du[5]
        D(Ca4CaMB_sl) ~ 1e-3*(rcn24B_sl + rcn4B_sl),                # du[6]
        # CaMKII equations
        D(Pb2_sl) ~ 1e-3*(rcnCKib2_sl - rcnCKb2b_sl + rcnCKt2b2_sl),     # du[7]
        D(Pb_sl) ~ 1e-3*(rcnCKib_sl + rcnCKb2b_sl - rcnCKbt_sl),         # du[8]
        D(Pt_sl) ~ 1e-3*(rcnCKbt_sl-rcnCKta_sl-rcnCKtt2_sl),             # du[9]
        D(Pt2_sl) ~ 1e-3*(rcnCKtt2_sl-rcnCKt2a_sl-rcnCKt2b2_sl),         # du[10]
        D(Pa_sl) ~ 1e-3*(rcnCKta_sl+rcnCKt2a_sl-rcnCKai_sl),             # du[11]
        # CaN equations
        D(Ca4CaN_sl) ~ 1e-3*(rcnCa4CaN_sl - rcn0CaN_sl - rcn2CaN_sl - rcn4CaN_sl),      # du[12]
        D(CaMCa4CaN_sl) ~ 1e-3*(rcn0CaN_sl - rcn02CaN_sl),                        # du[13]
        D(Ca2CaMCa4CaN_sl) ~ 1e-3*(rcn2CaN_sl+rcn02CaN_sl-rcn24CaN_sl),              # du[14]
        D(Ca4CaMCa4CaN_sl) ~ 1e-3*(rcn4CaN_sl+rcn24CaN_sl),                       # du[15]
        ## For adjusting Ca buffering in EC coupling model
        JCaSL ~ 1e-3*(2*CaMKIItot*(rcnCKtt2_sl-rcnCKb2b_sl) - 2*(rcn02_sl+rcn24_sl+rcn02B_sl+rcn24B_sl+rcnCa4CaN_sl+rcn02CaN_sl+rcn24CaN_sl))   # [uM/msec]
    ]
    
    eqs = append!(fluxes_eqs, CaMSL_eqs)

    @named sys = ODESystem(eqs, t)

    return sys
end

