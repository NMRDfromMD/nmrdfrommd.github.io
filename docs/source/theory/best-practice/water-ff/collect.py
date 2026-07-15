#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt


# SPCE

data_path = "/media/simon/data1/data-nmrdfrommd/SPCE/analysis/results/"

co = 14
N = 4000
T = [280, 290, 300, 310, 320]

out = []

for T0 in T:
    folder = f"N{N}co{co}T{T0}/"
    inter = data_path + f"N{N}_cutoff{co}_T{T0}_inter_molecule.npy" 
    intra = data_path + f"N{N}_cutoff{co}_T{T0}_intra_molecule.npy" 
    dataT = np.load(inter, allow_pickle=True).item()
    dataR = np.load(intra, allow_pickle=True).item()
    n_intra = dataT["n_intra"]
    n_inter = dataT["n_inter"]

    f = np.real(dataT["f"])
    R1T = np.real(dataT["R1"])
    R1T_std = np.real(dataT.get("R1_std"))
    R1T_error = R1T_std/np.sqrt(n_inter)

    R1R = np.real(dataR["R1"])
    R1R_std = np.real(dataR.get("R1_std"))
    R1R_error = R1R_std/np.sqrt(n_intra)

    R1 = R1T[0] + R1R[0]
    R1_error = np.sqrt(R1T_error[0]**2 + R1R_error[0]**2)

    T1 = 1/R1
    T1_error = R1_error / R1**2

    out.append([T0, T1, T1_error, R1T[0], R1T_error[0], R1R[0], R1R_error[0]])

    if T0 == 300:

        jump = 10

        t = np.real(dataT["t"])[::jump]
        GT = np.real(dataT["C"][0][::jump])
        GR = np.real(dataR["C"][0][::jump])
        GT_error = np.real(dataT["gij_std"][0][::jump])/np.sqrt(n_inter)
        GR_error = np.real(dataR["gij_std"][0][::jump])/np.sqrt(n_intra)
        
        np.savetxt("GR_T300K_SPCE.dat", np.vstack([t, GR, GR_error]).T, header="t Gij Gij_error", fmt="%.8e")
        np.savetxt("GT_T300K_SPCE.dat", np.vstack([t, GT, GT_error]).T, header="t Gij Gij_error", fmt="%.8e")

out = np.array(out)

np.savetxt("T1_SPCE.dat", out, header="T T1 T1_error R1T R1T_error R1T R1R_error", fmt="%.8e")
