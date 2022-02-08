import pandas as pd
import numpy as np

gas_density=10.0**3

n_particles=250
avs=np.logspace(-7,4,n_particles)

av_per_column=6.289E-22
column_densities=avs/av_per_column

x_values=column_densities/gas_density

gas_temp=50.0
dust_temp=20.0
chi=0.0
particle_type="P"


columns= ["Particle #","x (cm)","y (cm)","z (cm)","n_H (c^m-3)",
        "T_gas (K)","T_dust (K)"," chi (Draine)"," Type (I|P|D)"]


cloud=pd.DataFrame({
                columns[0]:range(0,n_particles),
                columns[1]:x_values,
                })

cloud[columns[2]]=0.0
cloud[columns[3]]=0.0

cloud[columns[4]]=gas_density

cloud[columns[5]]=gas_temp
cloud[columns[6]]=dust_temp
cloud[columns[7]]=chi
cloud[columns[8]]=particle_type
cloud.to_csv("Clouds/1e3.cloud",index=False)