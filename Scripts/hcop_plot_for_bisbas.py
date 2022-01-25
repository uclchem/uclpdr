import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use("thesis")
df=pd.read_csv("Outputs/var_Avenh_490.csv")

fig,ax=plt.subplots(figsize=(4.5,4.5),tight_layout=True)
ax.plot(df["n_H"],df["HCO+_abun"])
ax.set(xscale="log",yscale="log",
       xlabel="Density / cm$^{-3}$",
       ylabel="HCO$^+$")
fig.savefig("hcop_dens.png",dpi=300)

fig,ax=plt.subplots(figsize=(4.5,4.5),tight_layout=True)
ax.plot(df["Av"],df["HCO+_abun"])
ax.set(xscale="log",yscale="log",
       xlabel="Av / Mag",
       ylabel="HCO$^+$")
fig.savefig("hcop_av.png",dpi=300)