from plotting_functions import *
import matplotlib.pyplot as plt
######################################################################
#Settings
#Set variables to determine what we'll plot and how to arrange it
#####################################################################
file_prefix="ngc253_xrays" #prefix (including path) to output of model
ray=0 #pick a ray remembering to count from 0

dataframe_file="ngc253_xray_tests/ngc253_gmc5_full-xdr.csv" #filename to store uclpdr output in a machine readable form
create_dataframe=False #flag for whether dataframe exists

output_file="ngc253_xray_tests/ngc253_gmc5_full-xdr_heating.pdf"#dataframe_file.replace(".csv",".pdf")

figure_shape=(2,2) #nrows,ncols for grid of subplots
figure_size=(16,9) #size in inches, just sets ratio for pdf outputs

#fig columns should be a list of lists or strings.
#where you put a string, any column with that string in name will be plotted
#where you put a list, those specific columns are plotted
#_emis,_abun,_cool,_heat,_tau have been appended to columns so it's easy to get all of one type
# fig_columns=[
# 	"emis",["T_g","T_d"],["H+_abun","C_abun","C+_abun","H2_abun","H_abun"],"cool"
# 	]ls 

fig_columns=[
	["T_g","T_d"],["H+_abun","C_abun","C+_abun","H2_abun","H_abun","CO_abun"],"heat","cool"
	]


#these functions are in plotting_functions and will rename a given type of column
renamers=[
	self,abundance_rename,self,self
	]

#axis settings can be put in dictionaries so make a list of them to set random properties of each axis
axis_settings=[
				{"ylabel":"Temperature / K"},
				{"ylabel":"Fractional Abundance"},
                {"ylabel":"$\Lambda$ / erg cm$^{-3}$ s$^{-1}$","ylim":(1e-30,1e-16)},
				{"ylabel":"$\Lambda$ / erg cm$^{-3}$ s$^{-1}$","ylim":(1e-30,1e-16)}]


#save some effort by having an additional dictionary for settings you want for all axes
all_axis_settings={"xlabel":"A$_V$","xscale":"log","yscale":"log"}
########################################################################
# Load output
# Here we load the uclpdr raw output, save our version or load our version
# if we've run this script before.
########################################################################

#first load or create output file dataframe
if create_dataframe:
	#this function loads all uclpdr outputs into one table
	output_df=create_output_dataframe(file_prefix,ray)
	clear_raw_files(file_prefix)
	#we can write that table out and then we don't need the original .out files
	output_df.to_csv(dataframe_file,index=False)
else:
	output_df=pd.read_csv(dataframe_file)


##########################################################################
# Plotting
# Use output_df and setting to create and format a figure
##########################################################################

fig,axes=plt.subplots(figure_shape[0],figure_shape[1],figsize=figure_size,tight_layout=True)
axes=axes.flatten() #give a 1d array of axes instead of nrows.ncols array

#now loop through column choices and make plots
for i,cols in enumerate(fig_columns):
	renamer=renamers[i]
	print(cols)
	if type(cols)==str:
		#this gives a list of columns in table that contain the cols string
		cols=[x for x in output_df.columns if cols in x]
	for col in cols:
		axes[i].plot(output_df["Av"],output_df[col],label=renamer(col))


	axes[i].legend()
	axes[i].set(**axis_settings[i],**all_axis_settings)

fig.savefig(output_file)