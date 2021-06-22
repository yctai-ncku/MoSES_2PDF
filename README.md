**************************************************************************
##		*  *  *  *  *      MoSES_2PDF     *  *  *  *  *		
**************************************************************************
#### There are two modes:

- "adNOC_2Pb_landslide" is for the mode of landslide type.

- "adNOC_2Pb_inflow"    is for the mode of inflow type.

## Requirements

MoSES_2PDF was developed and tested on Ubuntu 18.04 LTS, in which the CUDA Version 11.0 is installed. 
The execution has been tested under CUDA Version 11.0 and 11.2.

- ##### PS How to install CUDA on Linux
```
sudo apt install nvidia-cuda-toolkit
```

## First step: Compilation

Use the makefile provided in the package 

(Remark: It works only with GCC-compiler. The IBM and PGI compilers do not work in the current version.)

- #### for MoSES_2PDF landslide type:
```
cd ./MoSES_2PDF_landslide
```
```
make 
```
- #### for MoSES_2PDF inflow type:
```
cd ./MoSES_2PDF_inflow
```

```
make 
```

## Second step: Build the result directory

- #### for MoSES_2PDF landslide type:
```
	mkdir ./MoSES_2PDF_landslide/result2Pb
	mkdir ./MoSES_2PDF_landslide/result2Pb/ascFile
	mkdir ./MoSES_2PDF_landslide/result2Pb/openGLFile
```
- #### for MoSES_2PDF inflow type:
```
	mkdir ./MoSES_2PDF_inflow/result2Pb
	mkdir ./MoSES_2PDF_inflow/result2Pb/ascFile
	mkdir ./MoSES_2PDF_inflow/result2Pb/openGLFile
```

## Third step: Parameter setting

All the parameters (simulation duration, material parameters and locations of the topography data) can be found and specified in "par_List".

- #### for MoSES_2PDF landslide type:
   in ./MoSES_2PDF_landslide/par_list

- #### for MoSES_2PDF inflow type:
    in ./MoSES_2PDF_inflow/par_list

## Fourth step: Execution

- #### execution for MoSES_2PDF landslide type:
```
    cd ./MoSES_2PDF_landslide
    ./MoSES_2PDF_landslide
```
- #### execution for MoSES_2PDF inflow type:
```
  cd ./MoSES_2PDF_inflow
    ./MoSES_2PDF_inflow
```
## TOPOGRAPHY DATA: (in directory "Data")

Following the government regulations, only the DEMs with a resolution of 20 m are provided in this package for test.


## DataProcess

- #### With plotting tool written in phthon
In directory "plotResult": The python code "plot_2PbGPU.py" can plot the computed results in directory "result2Pb".

- #### For GIS (asc-format)
In directory "toGISFile": The code "toGIS" can transform the computed results in directory "result2Pb" to directory "result2Pb/ascFile" for illustration in QGIS or ArcGIS.

- #### Illustration in ANSI-Platform
In directory "toANSIPFile": The code "toANSIP" can transform the computed results in directory "result2Pb" to directory "result2Pb/openGLFile" for illustration in ANSI-platform.

## Demo video

- #### MoSES_2PDF_demo:
https://youtu.be/9E9veNhrSME

- #### ANSIP_demo:
https://youtu.be/m0tNiK7Di6U

## Help

- In case of any question by applying the above code(s), please feel free to contact Dr. Yih-Chin Tai (email: yctai@ncku.edu.tw).

- And any suggestion(s) or collaboration for extending the current code is welcome.

Prof. Dr.-Ing Yih-Chin Tai <br>
Dept. Hydraulic and Ocean Engineering<br>
National Cheng Kung University, Tainan, Taiwan<br>
First version:  2021/04/07 (modified on 2021/06/22)
