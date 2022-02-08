from plotting_functions import *
from sys import argv

######################################################################
#Settings
#Set variables to determine what we'll plot and how to arrange it
#####################################################################
file_prefix=argv[1] #prefix (including path) to output of model


ray=0 #pick a ray remembering to count from 0

dataframe_file="Meijerink/"+file_prefix+".csv" #filename to store uclpdr output in a machine readable form


#this function loads all uclpdr outputs into one table
output_df=create_output_dataframe(file_prefix,ray)

#we can write that table out and then we don't need the original .out files
output_df.to_csv(dataframe_file,index=False)

#delete the files UCL-PDR produces now you have the full output dataframe
clear_raw_files(file_prefix)