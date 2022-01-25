import pandas as pd
import numpy as np
x,y,z,dens=np.loadtxt("var_Avenh_490.dat",
                         delimiter=",",unpack=True)
n_particles=len(x)

gas_temp=50.0
dust_temp=20.0
chi=0.0
particle_type="P"

columns= ["Particle #","x (cm)","y (cm)","z (cm)","n_H (c^m-3)",
        "T_gas (K)","T_dust (K)"," chi (Draine)"," Type (I|P|D)"]


cloud=pd.DataFrame({
                columns[0]:range(0,n_particles),
                columns[1]:x*3.086e+18,
                })

cloud[columns[2]]=0.0
cloud[columns[3]]=0.0

cloud[columns[4]]=dens

cloud[columns[5]]=gas_temp
cloud[columns[6]]=dust_temp
cloud[columns[7]]=chi
cloud[columns[8]]=particle_type
cloud.to_csv("Clouds/var_Avenh_490.cloud",index=False)