import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

N = 12
T = 100.0
B = 0.01

FileName = "./../data/specFn_N" + str(N) + "_T" + "{:.6f}".format(T) + "_B" + "{:.6f}".format(B) + ".csv"

df = pd.read_csv(FileName)
frequencies = df['frequency'].values
spectralFunction = df['spectralFunction'].values

fig = plt.figure()
plt.plot(frequencies, spectralFunction)
plt.xscale("log")
plt.yscale("log")
plt.show()