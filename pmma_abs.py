import matplotlib.pyplot as plt


abs = [10.55, 18.23, 24.6, 36.07, 39.7, 42.6, 43.69, 45.77, 52.97, 61.48, 66.44, 70.39, 79.11]
wl = [700.,  600.,  550.,  500.,  450.,  400.,  390.,  380.,  370.,  350.,  320., 310.,  300.]

plt.plot(wl,abs,"-o")

plt.title("PMMA absorption length")
plt.xlabel("wavelength, nm")
plt.ylabel("Mean free path, mm")

plt.show()