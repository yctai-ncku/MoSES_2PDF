**************************************************************************
		* * * * *     MoSES_2PDF (2 modes)     * * * * *		
**************************************************************************
"adNOC_2Pb_landslide" is for the mode of landslide type.

"adNOC_2Pb_inflow"    is for the mode of inflow type.

## Requirements

MoSES_2PDF was developed and tested on Ubuntu 18.04 LTS, in which the CUDA Version 11.2 is installed.

#### PS How to install CUDA on Linux
```
sudo apt install nvidia-cuda-toolkit
```

## First step: Compile to run  

Use the makefile provided in the package to compile

for MoSES_2PDF landslide type:
```
cd ./MoSES_2PDF_landslide
```
```
make 
```
for MoSES_2PDF inflow type:
```
cd ./MoSES_2PDF_inflow
```

```
make 
```

## Second step: build the result directory
	mkdir ./MoSES_2PDF_landslide/result2Pb
	mkdir ./MoSES_2PDF_landslide/result2Pb/ascFile
	mkdir ./MoSES_2PDF_landslide/result2Pb/openGLFile
	
	mkdir ./MoSES_2PDF_inflow/result2Pb
	mkdir ./MoSES_2PDF_inflow/result2Pb/ascFile
	mkdir ./MoSES_2PDF_inflow/result2Pb/openGLFile

## Third step (Mode-I): run the code for flows in landslide type:
    cd ./MoSES_2PDF_landslide
    ./MoSES_2PDF_landslide
## Third step (Mode-II): run the code for flows in inflow type:
    cd ./MoSES_2PDF_inflow/
    ./MoSES_2PDF_inflow


## PARAMETER SETTING: 
par_List


## TOPOGRAPHY DATA: (in directory "Data")

Following the government regulations, only the DEMs with a resolution of 20 m are provided.


## DataProcess

(1) In directory "plotResult": The python code "plot_2PbGPU.py" can plot the computed results in directory "result2Pb".

(2) In directory "toGISFile": The code "toGIS" can transform the computed results in directory "result2Pb" to directory "result2Pb/ascFile" for illustration in QGIS or ArcGIS.

(3) In directory "toANSIPFile": The code "toANSIP" can transform the computed results in directory "result2Pb" to directory "result2Pb/openGLFile" for illustration in ANSI-platform.

## Demo video

MoSES_2PDF_demo:
https://youtu.be/9E9veNhrSME

ANSIP_demo:
https://youtu.be/m0tNiK7Di6U

## Help

In case of any question by applying the above code(s), please feel free to contact Dr. Yih-Chin Tai.

And any suggestion(s) or collaboration for extending the current code is welcome.

Prof. Dr.-Ing Yih-Chin Tai (email: yctai@ncku.edu.tw)

Dept. Hydraulic and Ocean Engineering

National Cheng Kung University, Tainan, Taiwan

First version:  2021/04/07 (modified on 2021/06/21)
