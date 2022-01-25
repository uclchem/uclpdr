from plotting_functions import *
import matplotlib.pyplot as plt
import seaborn as sns
######################################################################
#Settings
#Set variables to determine what we'll plot and how to arrange it
#####################################################################
file_prefix="fixed_cooling" #prefix (including path) to output of model
ray=0 #pick a ray remembering to count from 0

dataframe_file="Benchmarks/bisbas_output.csv" #filename to store uclpdr output in a machine readable form
create_dataframe=False #flag for whether dataframe exists


figure_shape=(2,2) #nrows,ncols for grid of subplots
figure_size=(16,9) #size in inches, just sets ratio for pdf outputs

#fig columns should be a list of lists or strings.
#where you put a string, any column with that string in name will be plotted
#where you put a list, those specific columns are plotted
#_emis,_abun,_cool,_heat,_tau have been appended to columns so it's easy to get all of one type
# fig_columns=[
# 	"emis",["T_g","T_d"],["H+_abun","C_abun","C+_abun","H2_abun","H_abun"],"cool"
# 	]

fig_columns=[
	["C+_abun","C_abun","CO_abun"],["C+_abun","C_abun","CO_abun"],["H+_abun","H2_abun","H_abun"],["T_g","T_d"]
	]


#these functions are in plotting_functions and will rename a given type of column
renamers=[
	emissivity_rename,self,abundance_rename,self
	]

#axis settings can be put in dictionaries so make a list of them to set random properties of each axis
axis_settings=[{"ylabel":"Fractional Abundance"},
				{"ylabel":"Fractional Abundance","xlim":(1e3,1e6),"ylim":(1e-8,3e-4)},
				{"ylabel":"Fractional Abundance"},
				{"ylabel":"Temperature / K","xlim":(1e-3,1e3)}]


#save some effort by having an additional dictionary for settings you want for all axes
all_axis_settings={"xlabel":"n_H","xscale":"log","yscale":"log"}
########################################################################
# Load output
# Here we load the uclpdr raw output, save our version or load our version
# if we've run this script before.
########################################################################

dataframe_file="Benchmarks/bisbas_output.csv" #filename to store uclpdr output in a machine readable form
output_df=pd.read_csv(dataframe_file)
output_df["Mg"]="With MG"

dataframe_file="Benchmarks/bisbas_no_mg_output.csv"
output_df=output_df.append(pd.read_csv(dataframe_file))
output_df["Mg"]=output_df["Mg"].fillna("No Mg")


fig,axes=plt.subplots(figure_shape[0],figure_shape[1],figsize=figure_size,tight_layout=True)
axes=axes.flatten() #give a 1d array of axes instead of nrows.ncols array

#now loop through column choices and make plots
for i,cols in enumerate(fig_columns):
	plot_df=output_df[["n_H","Mg"]+cols]
	plot_df=plot_df.melt(id_vars=["n_H","Mg"])
	sns.lineplot(data=plot_df,x="n_H",y="value",hue="variable",style="Mg",ax=axes[i])
	axes[i].legend()
	axes[i].set(**axis_settings[i],**all_axis_settings)
# ##########################################################################
# # Plotting
# # Use output_df and setting to create and format a figure
# ##########################################################################


# 	renamer=renamers[i]
# 	colors=sns.color_palette()

# 	if type(cols)==str:
# 		#this gives a list of columns in table that contain the cols string
# 		cols=[x for x in output_df.columns if cols in x]
# 	for col in cols:
# 		axes[i].plot(output_df["n_H"],output_df[col],label=renamer(col),ls="--",color=colours[i])


# 	axes[i].legend()
# 	axes[i].set(**axis_settings[i],**all_axis_settings)

# dataframe_file="Benchmarks/bisbas_no_mg_output.csv"
# output_df=pd.read_csv(dataframe_file)#filename to store uclpdr output in a machine readable form
# #now loop through column choices and make plots
# for i,cols in enumerate(fig_columns):
# 	renamer=renamers[i]
# 	colors=sns.color_palette()
# 	if type(cols)==str:
# 		#this gives a list of columns in table that contain the cols string
# 		cols=[x for x in output_df.columns if cols in x]
# 	for j,col in enumerate(cols):
# 		axes[i].plot(output_df["n_H"],output_df[col],label=renamer(col),color=colours[i])


	

fig.savefig("bisbas.pdf",type="PDF")