import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as mpl


schreiben = 0
## plot result type
## 1 : depth, 2 : phiS0, 3 : velocity Solid, 4 : velocity Fluid  
OutType = 1

## input Figure step
FigNr = 0
iistep = 60
iiend  = 62

dir_loc = '../../result2Pb/'
dir_out = '../../result2Pb/'

info_file = dir_loc+'Info.dat'
topo_file = dir_loc+'DEM.dat'
time_file = dir_loc+'Time.dat'

dataInfo = np.loadtxt(info_file,skiprows=1)
dataTopo = np.loadtxt(topo_file,skiprows=1)
dataTime = np.loadtxt(time_file)

## Information of time scale
TimeScle = math.sqrt(0.1/9.8)
timeStep = np.round(dataTime*TimeScle,2)

## Information of topography
NX = int(dataInfo[0])
NY = int(dataInfo[1])
dx = dataInfo[2]
dy = dataInfo[3]

TotalStep = dataInfo[6]

arraysize = int(NX*NY)

## Topography Data from DEM.plt
x_ini  = dataTopo[0:arraysize,0]
y_ini  = dataTopo[0:arraysize,1]
z_ini  = dataTopo[0:arraysize,2]
c_ini  = dataTopo[0:arraysize,3]
S1_ini = dataTopo[0:arraysize,4]
S2_ini = dataTopo[0:arraysize,5]


## plot topo information 
x_loc = x_ini.reshape([NX,NY])
y_loc = y_ini.reshape([NX,NY])
topo  = z_ini.reshape([NX,NY])
cvalue= c_ini.reshape([NX,NY])

flowIndex = np.zeros([371,222])


for ii in range(1,iiend,iistep):
    phys_input = "%03d" % ii+'.dat'
    phys_file  = dir_loc+phys_input
    dataphys   = np.loadtxt(phys_file,skiprows=1)
    
    H_ini  = dataphys[0:arraysize,0]
    Phi_ini= dataphys[0:arraysize,1]
    Us_ini = dataphys[0:arraysize,2]
    Uf_ini = dataphys[0:arraysize,3]
    Vs_ini = dataphys[0:arraysize,4]
    Vf_ini = dataphys[0:arraysize,5]
    
    H   =   H_ini.reshape([NX,NY])
    Phi = Phi_ini.reshape([NX,NY])
    Us  =  Us_ini.reshape([NX,NY])
    Uf  =  Uf_ini.reshape([NX,NY])
    Vs  =  Vs_ini.reshape([NX,NY])
    Vf  =  Vf_ini.reshape([NX,NY])
  
    h_indx = np.ones(H.shape,dtype=np.float)
    h_indx[(H[:,:]<0.01).nonzero()] = np.nan
    
    Mh       = (H/cvalue)*h_indx
    MPhi     = (Phi)*h_indx
    MspeedS = (np.sqrt(Us**2+Vs**2))*h_indx
    MspeedF = (np.sqrt(Uf**2+Vf**2))*h_indx
    
    flowIndex[np.where(H[:,:]>0.01)] = 1
    
    ## calculate colorbar range
    if ii==1:
        if OutType == 1:
            maxH = np.max(H/cvalue)
            minH = np.min(H/cvalue)
            if maxH%10!=0:
                maxH = int(maxH/10)*10+10
            if minH%10!=0:
                minH = int(minH/10)*10-10
        elif OutType == 2:
            maxH = 0.7#np.max(Phi)
            minH = 0#np.min(Phi)
        elif OutType == 3:
            maxH = 20 #np.max(np.sqrt(Us**2+Vs**2))
            minH = 0  #np.min(np.sqrt(Us**2+Vs**2))
        elif OutType == 4:
            maxH = 20 #np.max(np.sqrt(Uf**2+Vf**2))
            minH = 0 #np.min(np.sqrt(Uf**2+Vf**2))
        
        numTmp = (maxH-minH)/5
    
    ## topography for contour 
    fig = plt.figure(figsize=(NX/20,NY/20))
    ax = plt.gca()
    ax.set_aspect(1)
    TopoImage = plt.contour(x_loc,y_loc,topo, 100,colors='black')
    plt.xticks(())
    plt.yticks(())
    
    color_map = mpl.cm.rainbow

    
    if OutType == 1:
        Ch2 = plt.pcolormesh(x_loc, y_loc, Mh, shading='auto')
        figTXT = 'Fig_H'+str(ii-1)
    elif OutType == 2:
        Ch2 = plt.pcolormesh(x_loc, y_loc, MPhi, shading='auto')
        figTXT = 'Fig_PhiS0'+str(ii-1)
    elif OutType == 3:
        Ch2 = plt.pcolormesh(x_loc, y_loc, MspeedS, shading='auto')
        figTXT = 'Fig_SpeedS'+str(ii-1)
    elif OutType == 4:
        Ch2 = plt.pcolormesh(x_loc, y_loc, MspeedF, shading='auto')
        figTXT = 'Fig_SpeedF'+str(ii-1)
    
    cbar = plt.colorbar(shrink = 0.82,aspect = 20)
    cbar.ax.tick_params(labelsize=24)

    plt.text(0,max(y_ini)*1.02,'t = '+ str(timeStep[ii-1]) +' sec', fontdict={'size': 16, 'color': 'black'},backgroundcolor = 'w')
    plt.show()
    
    if schreiben == 1: 
        plt.savefig(dir_out+figTXT+'.png') 
    



