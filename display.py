import sys, os, csv

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

print("Running")
os.system("a.exe") #Run the programme
print("Simulation Finished")

data_reader = pd.read_csv('numerov.csv', header = None)

y = data_reader.iloc[0].values

x = np.linspace(-2, 2, 401)

plt.plot(x,y)
plt.show()
