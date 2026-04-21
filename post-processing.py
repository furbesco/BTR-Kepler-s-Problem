import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('kepler_output.csv')

t   = df["t"].to_numpy()
E   = df["l"].to_numpy() 
u   = df["u"].to_numpy() 
R   = df["R"].to_numpy()
v   = df["v"].to_numpy() 
phi   = df["phi"].to_numpy() 
x = df["x"].to_numpy()
y = df["y"].to_numpy()


plt.figure() 
plt.plot(x, y, label="Orbit")
plt.scatter([0], [0], s=20, label="Centre Mass") 
plt.gca().set_aspect("equal", adjustable="box")
plt.xlabel("x") 
plt.ylabel("y")
plt.title("1PN Orbit")
plt.legend() 
plt.tight_layout()
plt.show() 
