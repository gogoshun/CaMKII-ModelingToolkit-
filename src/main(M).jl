#=
    This file integrate all ODE files.
=#

include("camdyad_ODEfile(M).jl")
include("camsl_ODEfile(M).jl")
include("camcyt_ODEfile(M).jl")
include("camkii_ODEfile(M).jl")
include("ecc_ODEfile(M).jl")
include("bar_ODEfile(M).jl")

using ModelingToolkit
using DifferentialEquations
using Plots


sys_camdyad = get_camdyad_equations()
sys_camsl = get_camsl_equations()
sys_camcyt = get_camcyt_equations()
sys_camkii = get_camkii_equations()
sys_ecc =  get_ecc_equations()
sys_bar = get_bar_equations()

sys1 = extend(sys_camkii, sys_camsl)
sys2 = extend(sys1, sys_camcyt)
sys3 = extend(sys2, sys_camdyad)
sys4 = extend(sys3, sys_ecc)
sys = extend(sys4, sys_bar)

# sys = structural_simplify(eqs)


@unpack Na_m, Na_h, Na_j, ICa_HH4, ICa_HH5, ICa_HH6, ICa_HH7, Itos_x, Itos_y, Itof_x, Itof_y, Ikr, IKs, RyR_R, RyR_O, RyR_I, NaBj, NaBsl,  
        TnCL, TnCHc, TnCHm, CaM, Myosin_ca, Myosin_mg, SRB, SLLj, SLLsl, SLHj, SLHsl, Csqn, Ca_sr, Naj, Nasl, Nai, Ki, Ca_j, Ca_sl, Cai, Vm, 
        Itos_r, influx_LTCC, influx_PMCA, influx_NCX, influx_ICa, Na_late_h, CNa2, CNa1, ONa, IFNa, I1Na, CNa3, ICNa2, ICNa3, LONa, LCNa1, LCNa2, LCNa3, 
        C2_m1j, C1_m1j, I1Ca_m1j, I2Ca_m1j, I1Ba_m1j, I2Ba_m1j, C2_m2j, C1_m2j, I1Ca_m2j, I2Ca_m2j, I1Ba_m2j, I2Ba_m2j, C2_m1sl, C1_m1sl, I1Ca_m1sl, 
        I2Ca_m1sl, I1Ba_m1sl, I2Ba_m1sl, C2_m2sl, C1_m2sl, I1Ca_m2sl, I2Ca_m2sl, I1Ba_m2sl, I2Ba_m2sl, IKs_x, IKs1_y, Iss, IKs2_y,  # ecc_ODEfile
        CaM_dyad, Ca2CaM_dyad, Ca4CaM_dyad, CaMB_dyad, Ca2CaMB_dyad, Ca4CaMB_dyad, Pb2_dyad, Pb_dyad, 
        Pt_dyad, Pt2_dyad, Pa_dyad, Ca4CaN_dyad, CaMCa4CaN_dyad, Ca2CaMCa4CaN_dyad, Ca4CaMCa4CaN_dyad,                              # camdyad_ODEfile
        CaM_sl, Ca2CaM_sl, Ca4CaM_sl, CaMB_sl, Ca2CaMB_sl, Ca4CaMB_sl, Pb2_sl, Pb_sl, 
        Pt_sl, Pt2_sl, Pa_sl, Ca4CaN_sl, CaMCa4CaN_sl, Ca2CaMCa4CaN_sl, Ca4CaMCa4CaN_sl,                                            # camsl_ODEfile
        CaM_cyt, Ca2CaM_cyt, Ca4CaM_cyt, CaMB_cyt, Ca2CaMB_cyt, Ca4CaMB_cyt, Pb2_cyt, Pb_cyt, 
        Pt_cyt, Pt2_cyt, Pa_cyt, Ca4CaN_cyt, CaMCa4CaN_cyt, Ca2CaMCa4CaN_cyt, Ca4CaMCa4CaN_cyt,                                     # camcyt_ODEfile
        LCC_PKAp, LCC_CKdyadp, RyR2809p, RyR2815p, PLBT17p, LCC_CKslp,                                                              # camkii_ODEfile
        LR, LRG, RG, b1AR_S464, b1AR_S301, GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP, PDEp, cAMPtot, RC_I, RCcAMP_I, 
        RCcAMPcAMP_I, RcAMPcAMP_I, PKACI, PKACI_PKI, RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII,                         # bar_ODEfile
        PKACII_PKI, I1p_PP1, I1ptot, PLBp, PLMp, LCCap, LCCbp, RyRp, TnIp, KS79, KS80, KSp, CFTRp, KURp = sys


oprob = ODEProblem(sys, 
                    [Na_m => 1.94e-3, Na_h => 0.981, Na_j => 0.987, ICa_HH4 => 7.02e-6, ICa_HH5 => 1.00068, ICa_HH6 => 2.7e-2, ICa_HH7 => 1.6e-2, Itos_x => 2.02e-3, Itos_y => 0.99, Itof_x => 2.02e-3, 
                    Itof_y => 0.9992, Ikr => 1.11e-2, IKs => 7.37e-3, RyR_R => 0.698, RyR_O => 4.24e-6, RyR_I => 1.84e-6, NaBj => 3.993, NaBsl => 0.87, TnCL => 9.26e-3, TnCHc => 0.118, 
                    TnCHm => 1.03e-2, CaM => 2.53e-4, Myosin_ca => 1.989e-3, Myosin_mg => 0.138, SRB => 2.26e-3, SLLj => 2.2e-2, SLLsl => 1.35e-2, SLHj => 0.127, SLHsl => 0.142, Csqn => 1.177, 
                    Ca_sr => 0.503, Naj => 11.182, Nasl => 11.182, Nai => 11.182, Ki => 134.99, Ca_j => 5.34e-4, Ca_sl => 1.46e-4, Cai => 9.12e-5, Vm => -83.632, Itos_r => 0.946, 
                    influx_LTCC => 5.59e4, influx_PMCA => -3.38e4, influx_NCX => -3.096e5, influx_ICa => 2.875e5, Na_late_h => 0.222, CNa2 => 0.105, CNa1 => 1.92e-3, ONa => 4.15e-5, IFNa => 0.303, I1Na => 0.566, 
                    CNa3 => 1.01e-2, ICNa2 => 7.01e-5, ICNa3 => 8.62e-8, LONa => 1.47e-4, LCNa1 => 2.64e-6, LCNa2 => 1.82e-8, LCNa3 => 2.17e-11, C2_m1j => 0.939, C1_m1j => 2.71e-5, I1Ca_m1j => 9.17e-5, 
                    I2Ca_m1j => 6.71e-4, I1Ba_m1j => 4.99e-5, I2Ba_m1j => 5.97e-2, C2_m2j => 0.939, C1_m2j => 2.71e-5, I1Ca_m2j => 9.13e-5, I2Ca_m2j => 6.69e-4, I1Ba_m2j => 4.99e-5, I2Ba_m2j => 5.97e-2, C2_m1sl => 0.94, 
                    C1_m1sl => 2.71e-5, I1Ca_m1sl => 3.88e-6, I2Ca_m1sl => 2.81e-5, I1Ba_m1sl => 5.00e-5, I2Ba_m1sl => 5.98e-2, C2_m2sl => 0.94, C1_m2sl => 2.71e-5, I1Ca_m2sl => 4.17e-6, I2Ca_m2sl => 3.02e-5, I1Ba_m2sl => 5.00e-5, 
                    I2Ba_m2sl => 5.98e-2, IKs_x => 7.37e-3, IKs1_y => 0.99, Iss => 7.37e-3, IKs2_y => 0.995, 
                    CaM_dyad => 388.68, Ca2CaM_dyad => 13.02, Ca4CaM_dyad => 8.68e-3, CaMB_dyad => 0.0, Ca2CaMB_dyad => 0.0, Ca4CaMB_dyad => 0.0, Pb2_dyad => 0.67, Pb_dyad => 7.09e-2, 
                    Pt_dyad => 2.42e-5, Pt2_dyad => 8.73e-9, Pa_dyad => 3.37e-9, Ca4CaN_dyad => 7.35e-5, CaMCa4CaN_dyad => 2.43e-3, Ca2CaMCa4CaN_dyad => 0.013, Ca4CaMCa4CaN_dyad => 3.6,
                    CaM_sl => 4.42e-2, Ca2CaM_sl => 7.34e-5, Ca4CaM_sl => 8.89e-9, CaMB_sl => 2.44, Ca2CaMB_sl => 11.86, Ca4CaMB_sl => 4.38e-4, Pb2_sl => 1.47e-5, Pb_sl => 6.31e-6, 
                    Pt_sl => 6.60e-8, Pt2_sl => 7.37e-13, Pa_sl => 4.37e-9, Ca4CaN_sl => 5.22e-4, CaMCa4CaN_sl => 1.98e-6, Ca2CaMCa4CaN_sl => 5.02e-6, Ca4CaMCa4CaN_sl => 1.43e-3,
                    CaM_cyt => 4.4e-2, Ca2CaM_cyt => 4.11e-5, Ca4CaM_cyt => 6.17e-10, CaMB_cyt => 4.179, Ca2CaMB_cyt => 1.11, Ca4CaMB_cyt => 1.61e-5, Pb2_cyt => 8.23e-6, Pb_cyt => 4.15e-8,
                    Pt_cyt => 2.29e-13, Pt2_cyt => 2.47e-18, Pa_cyt => 1.53e-14, Ca4CaN_cyt => 1.17e-4, CaMCa4CaN_cyt => 4.39e-7, Ca2CaMCa4CaN_cyt => 1.59e-7, Ca4CaMCa4CaN_cyt => 1.49e-6,
                    LCC_PKAp => 16.454, LCC_CKdyadp =>16.934 , RyR2809p => 297.36, RyR2815p => 76.985, PLBT17p => 0.614, LCC_CKslp => 8.66e-6,
                    LR => -6.7e-36, LRG => 2.46e-34, RG => 4.8e-4, b1AR_S464 => 5.97e-35, b1AR_S301 => 6.48e-4, GsaGTPtot => 9.6e-3, GsaGDP => 6.21e-4, Gsby => 0.01, AC_GsaGTP => 1.42e-3, PDEp => 2.22e-3, 
                    cAMPtot => 1.023, RC_I => 0.804, RCcAMP_I => 0.142, RCcAMPcAMP_I => 4.48e-3, RcAMPcAMP_I => 0.229, PKACI => 8.55e-2, PKACI_PKI => 0.144, RC_II => 0.051, RCcAMP_II => 8.99e-3, RCcAMPcAMP_II => 2.84e-4, 
                    RcAMPcAMP_II => 5.77e-2, PKACII => 2.15e-2, PKACII_PKI => 3.62e-2, I1p_PP1 => 7.27e-2, I1ptot => 7.28e-2, PLBp => 8.454, PLMp => 5.6, LCCap => 5.49e-3, LCCbp => 6.27e-3, RyRp => 2.76e-2, 
                    TnIp => 4.389, KS79 => 1.53e-3, KS80 => 1.53e-3, KSp => 1.84e-3, CFTRp => 4.06e-3, KURp => 1.09e-2], 
                   (0.0, 1e3))

sol = solve(oprob, saveat=0.01, abstol = 1e-9, reltol = 1e-9)

plot(sol, idxs=Cai, linewidth=3, title="Calcium Transient", xlabel="Time(s)", ylabel="[Ca2+](mM)")

plot(sol, idxs=Vm, linewidth=3, title="Action Potential", xlabel="Time(s)", ylabel="Voltage (mV)")
