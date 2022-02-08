from plotting_functions import *
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
plt.style.use("thesis")
######################################################################
#Settings
#Set variables to determine what we'll plot and how to arrange it
#####################################################################

#get GMC 5 abundances
abunds=pd.read_csv("/home/jon/Documents/zupcx4/CR_project/data/posterior_abundance_summaries.csv")
abunds=abunds.loc[2:2].melt()

abunds["variable"]=abunds["variable"].str.replace("<lambda_0>","")
abunds["variable"]=abunds["variable"].str.replace("<lambda_1>","")
abunds["variable"]=abunds["variable"].str.replace("median","")
abunds=abunds.groupby("variable").agg(["min","max"]).reset_index()

abunds.columns=["variable","min","max"]


models={"F$_x$=0.05":"./ngc253_xray_tests/full-network_gmc5_xrays_4e7.csv",
 "F$_x$=0.0":"./ngc253_xray_tests/full-network_gmc5_no_xrays.csv"}

styles=dict(zip(models.keys(),["-","--"]))
xmax=20

fig,ax=plt.subplots(1,1,figsize=(4,4),tight_layout=True)
colors=sns.color_palette("colorblind",n_colors=2)
for model,path in models.items():

    output_df=pd.read_csv(path)
    for i,species in enumerate(["H3O+","SO"]):
        label= "" if i==1 else model
        ax.plot(output_df["Av"],output_df[species+"_abun"]
                ,color=colors[i],ls=styles[model],label=label,zorder=3)
        spec,min_val,max_val=abunds.loc[abunds["variable"]==species].values[0]
        if model==list(models.keys())[0]:
                ax.axhspan(min_val,max_val,xmax=xmax,color="white",alpha=1.0,zorder=2-i)

                ax.axhspan(min_val,max_val,xmax=xmax,color=colors[i],alpha=.4,label=species,zorder=2-i)
ax.legend()
ax.set(yscale="log",xlim=(0,xmax),xlabel="A$_V$ / mag",ylabel="Fractional Abundance")
fig.savefig("/home/jon/Documents/zupcx4/CR_project/Paper/xray_model.pdf")
