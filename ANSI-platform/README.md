# ANSI-platform
Advanced Numerical Scenario Illustration Platform using OpenGL
## Requirements
* glew
* glfw3
* sdl2-image
### How to install on Linux
```
sudo apt-get install libglfw3-dev libglfw3
sudo apt-get install libglew-dev
sudo apt-get install libglm-dev
sudo apt-get install libsdl2-image-dev
pkg-config --libs --static glew
pkg-config --libs --static glfw3
```
## Compiling and running a particular demo
```
make all
```
## Compile
```
make compile
```
## Run
You can run the command `./output <folder-name>` where `<folder-name>` is the folder name of DEM files.

For example: 
```
./output terrain/D056_plt_run
```
## Control
### Light
```
'W': Move backward
'S': Move forward
'A': Move left
'D': Move right
```
### Terrain Position
```
'UP': Move upward
'DOWN': Move downward
'LEFT': Move left
'RIGHT': Move right
```
### Texture
```
'1': Grass
'2': Satellite
'3': Gray
```
### Transparency
```
']': 100%
'[': 50%
```
### Animation
```
'SPACE': Start
'BACKSPACE': Initialize
'X': Next frame
'Z': Prev frame
```
