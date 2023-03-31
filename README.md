# MeshGeodesy
Contains application and library for geodesic distance and path computation over triangle meshes. 
The code is mainly composed by:
* mgapp: Appication code that launches an XWindows to render mesh and its geodesic distance and path information through OpenGL.
* mglib: Library code that reads in a triangle mesh from file, and outputs geodesic distance and path of a selected vertex to all other vertices in the triangle mesh.

# Requirements
This project was tested under a Ubuntu system. Install the following dependency libraries: 
* Google glog
```sudo apt-get update```
```sudo apt-get install -y libgoogle-glog-dev```
* FLTK a cross platform GUI development library. Link to https://github.com/fltk/fltk/blob/master/README.txt for installation instructions.
