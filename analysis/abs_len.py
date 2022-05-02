import numpy as np
import matplotlib.pyplot as plt


pmma_wl            = np.array([700.,  600.,  550.,  500.,  450.,  400.,  390.,  380.,  370.,  350.,  320., 310.,  300. ]);
pmma_side_abslen   = np.array([10.55 , 18.23 , 24.6, 36.07, 39.7, 42.6, 43.69, 45.77, 52.97, 61.48, 66.44, 70.39, 79.11]);
pmma_bottom_abslen = np.array([ 0.01,  0.01,   0.01, 0.01, 1.31, 4.26, 14.98, 24.09, 28.56, 30.35, 32.39, 33.93, 37.32]);

opEn = np.array([
    2.06640, 2.10143, 2.13766, 2.17516, 2.21400, 2.25426, 2.29600, 2.33932, 2.38431, 2.43106,
    2.47968, 2.53029, 2.58300, 2.63796, 2.69531, 2.75520, 2.81782, 2.88335, 2.95200, 3.09960,
    3.54241, 4.13281 ])

AbsLen_PMMA = np.array([
    39.48, 48.25, 54.29, 57.91, 54.29, 33.40, 31.02, 43.43, 43.43, 41.36,
    39.48, 37.76, 36.19, 36.19, 33.40, 31.02, 28.95, 25.55, 24.13, 21.71,
    2.171, 0.434  ])

pmma_en = 1240./pmma_wl

# plt.plot(pmma_en, pmma_side_abslen, "-o",label="pmma_side")
# plt.plot(pmma_en, pmma_bottom_abslen, "-o",label="pmma_bottom")
plt.plot(opEn, AbsLen_PMMA, "-o")

plt.ylabel("Absorption length, m")
plt.xlabel("Photon energy, eV")
# plt.legend()
plt.show()