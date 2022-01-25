import pandas as pd
import numpy as np

n_particles=100
total_size=1.060051943E+18

gas_temp=50.0
dust_temp=20.0
chi=0.0
particle_type="P"


#define here a function that returns the density for any given x
def density_function(x):
    return 10.0**6.0

columns= ["Particle #","x (cm)","y (cm)","z (cm)","n_H (c^m-3)",
        "T_gas (K)","T_dust (K)"," chi (Draine)"," Type (I|P|D)"]

sizes=np.linspace(0,total_size,n_particles)

cloud=pd.DataFrame({
                columns[0]:range(0,n_particles),
                columns[1]:sizes,
                })

cloud[columns[2]]=0.0
cloud[columns[3]]=0.0

cloud[columns[4]]=cloud[columns[1]].map(density_function)

cloud[columns[5]]=gas_temp
cloud[columns[6]]=dust_temp
cloud[columns[7]]=chi
cloud[columns[8]]=particle_type
cloud.to_csv("Clouds/1e6.cloud",index=False)