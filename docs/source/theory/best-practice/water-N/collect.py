#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt



# TIP4P

data_path = "/media/simon/data1/data-nmrdfrommd/TIP4P/analysis/results/"

co = 14
N = [25, 38, 58, 88, 135, 207, 316, 482, 736, 1717, 2620, 4000, 6106]
T = 300

out = []

for N0 in N:

    all_R1T = []
    all_R1T_std = []
    all_R1T_error = []
    all_R1R = []
    all_R1R_std = []
    all_R1R_error = []

    all_GT = []
    all_GR = []
    all_GT_err = []
    all_GR_err = []

    nb_data_found = 0

    for idx in np.arange(1, 21):

        folder = f"N{N0}co{co}T{T}/"

        try:
            
            inter = data_path + f"N{N0}_cutoff{co}_T{T}_inter_molecule_idx{idx}.npy" 
            intra = data_path + f"N{N0}_cutoff{co}_T{T}_intra_molecule_idx{idx}.npy" 
            dataT = np.load(inter, allow_pickle=True).item()
            dataR = np.load(intra, allow_pickle=True).item()
            n_intra = dataT["n_intra"]
            n_inter = dataT["n_inter"]

            f = np.real(dataT["f"])
            all_R1T.append(np.real(dataT["R1"]))
            all_R1T_std.append(np.real(dataT.get("R1_std")))
            all_R1T_error.append(all_R1T_std[-1]/np.sqrt(n_inter))

            all_R1R.append(np.real(dataR["R1"]))
            all_R1R_std.append(np.real(dataR.get("R1_std")))
            all_R1R_error.append(all_R1R_std[-1]/np.sqrt(n_intra))

            t = np.real(dataT["t"])
            all_GT.append(np.real(dataT["C"]))
            all_GR.append(np.real(dataR["C"]))
            all_GT_err.append(np.real(dataT["gij_std"]))
            all_GR_err.append(np.real(dataR["gij_std"]))

            nb_data_found += 1

        except:
            
            pass
            # print("NOT FOUND", folder, idx)

    R1T = np.mean(all_R1T, axis = 0)
    R1R = np.mean(all_R1R, axis = 0)

    GT = np.mean(all_GT, axis = 0)
    GR = np.mean(all_GR, axis = 0)

    tau_T = np.trapz(GT[0], t) / GT[0][0]
    tau_R = np.trapz(GR[0], t) / GR[0][0]

    if nb_data_found == 1:
        R1T_std = np.mean(all_R1T_std, axis = 0)
        R1T_error = np.mean(all_R1T_error, axis = 0)
        R1R_std = np.mean(all_R1R_std, axis = 0)
        R1R_error = np.mean(all_R1R_error, axis = 0)
        GT_err = np.mean(all_GT_err, axis = 0)
        GR_err = np.mean(all_GR_err, axis = 0)
    else:
        R1T_std = np.std(all_R1T, axis = 0)
        R1R_std = np.std(all_R1R, axis = 0)
        R1T_error = R1T_std / np.sqrt(len(all_R1T))
        R1R_error = R1R_std / np.sqrt(len(all_R1R))
        GT_std = np.std(all_GT, axis = 0)
        GR_std = np.std(all_GR, axis = 0)
        GT_err = GT_std / np.sqrt(len(all_GT))
        GR_err = GR_std / np.sqrt(len(all_GR))

    R1 = R1T[0] + R1R[0]
    R1_error = np.sqrt(R1T_error[0]**2 + R1R_error[0]**2)

    T1 = 1/R1
    T1_error = R1_error / R1**2


    # ----- tau_T -----
    dt = np.diff(t)
    dt = np.append(dt, dt[-1])      # approximate spacing

    I_T = np.trapz(GT[0], t)
    sigma_I_T = np.sqrt(np.sum((GT_err[0] * dt)**2))

    tau_T = I_T / GT[0][0]

    tau_T_error = np.sqrt(
        (sigma_I_T / GT[0][0])**2 +
        (I_T * GT_err[0][0] / GT[0][0]**2)**2
    )

    # ----- tau_R -----
    I_R = np.trapz(GR[0], t)
    sigma_I_R = np.sqrt(np.sum((GR_err[0] * dt)**2))

    tau_R = I_R / GR[0][0]

    tau_R_error = np.sqrt(
        (sigma_I_R / GR[0][0])**2 +
        (I_R * GR_err[0][0] / GR[0][0]**2)**2
    )


    out.append([N0, R1, R1_error, R1T[0], R1T_error[0], R1R[0], R1R_error[0], tau_T, tau_R, tau_T_error, tau_R_error])

    # if T0 == 300:

    #     jump = 10

    #     t = np.real(dataT["t"])[::jump]
    #     GT = np.real(dataT["C"][0][::jump])
    #     GR = np.real(dataR["C"][0][::jump])
    #     GT_error = np.real(dataT["gij_std"][0][::jump])/np.sqrt(n_inter)
    #     GR_error = np.real(dataR["gij_std"][0][::jump])/np.sqrt(n_intra)
        
    #     np.savetxt("GR_T300K_TIP4P.dat", np.vstack([t, GR, GR_error]).T, header="t Gij Gij_error", fmt="%.8e")
    #     np.savetxt("GT_T300K_TIP4P.dat", np.vstack([t, GT, GT_error]).T, header="t Gij Gij_error", fmt="%.8e")

out = np.array(out)

np.savetxt("T1_TIP4P.dat", out, header="N0 R1 R1_error R1T R1T_error R1T R1R_error tauT tauR", fmt="%.8e")



