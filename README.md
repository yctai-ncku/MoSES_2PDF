**************************************************************************
		* * * * *     MoSES_2PDF mode     * * * * *		
**************************************************************************
"adNOC_2Pb_landslide" is the landslide type mode.
"adNOC_2Pb_inflow"    is the inflow type mode.

(PS If you need more cells in the mode, please feel free to contact Dr. Tai by email: yctai@ncku.edu.tw)

### Requirements
MoSES_2PDF was developed and tested on Ubuntu 18.04 LTS. It relies on CUDA Version 11.2.

### How to install on Linux
```
sudo apt install nvidia-cuda-toolkit
```
Use a makefile to compile

```
cd ./MoSES_2PDF_landslide
make 
```
or
```
cd ./MoSES_2PDF_inflow
make 
```

## build the result directory
	mkdir ./MoSES_2PDF_landslide/result2Pb
	mkdir ./MoSES_2PDF_landslide/result2Pb/ascFile
	mkdir ./MoSES_2PDF_landslide/result2Pb/openGLFile
	
	mkdir ./MoSES_2PDF_inflow/result2Pb
	mkdir ./MoSES_2PDF_inflow/result2Pb/ascFile
	mkdir ./MoSES_2PDF_inflow/result2Pb/openGLFile

## run MoSES_2PDF landslide type:
    cd ./MoSES_2PDF_landslide
    ./MoSES_2PDF_landslide
## MoSES_2PDF inflow type:
    cd ./MoSES_2PDF_inflow/
    ./MoSES_2PDF_inflow


## PARAMETER SETTING: par_List


## TOPOGRAPHY DATA: (in directory "Data")
Because of relevant regulations, We only provides terrain with a grid accuracy of 20 m.


## DataProcess
(1) "plotResult" has plot_2PbGPU.py. The python code can plot the "MoSES_2PDF" computation result.

(2) "toGIS" has ResulttoGIS. The code can transform the "MoSES_2PDF" computation result to QGIS or ArcGIS.

(3) "toANSIP" has ResulttoANSIP. The code can transform the "MoSES_2PDF" computation result to ANSI-platform.

## Demo video

MoSES_2PDF_demo:
https://youtu.be/9E9veNhrSME

ANSIP_demo:
https://youtu.be/m0tNiK7Di6U

## Help (email: yctai@ncku.edu.tw)
In case of any question by applying the above code(s), please feel free
to contact Dr. Yih-Chin Tai. And any suggestion(s) or collaboration for
extending the current code is welcome.

Associate Prof. Dr.-Ing Yih-Chin Tai
email: yctai@ncku.edu.tw
Dept. Hydraulic and Ocean Engineering
National Cheng Kung University, Taiwan

First version: 2021/04/07
