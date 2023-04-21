# Mesh Geodesy


Geodesic distance calculation is a powerful technique in mesh processing that allows measurement of the shortest path between two points on a 3D surface. Unlike Euclidean distance, which assumes a flat space, geodesic distance takes into account the curvature of the surface and provides a metric that is constrained to the geometry of the mesh. This makes it particularly useful for a variety of applications, in computer graphics, computational geometry, physical simulations, and other areas.

| Example | 
| ------- | 
| <img src="pictures/lucy3.gif" alt="alt text" height="500"> | 
| Progressive rendering of a 3D mesh with triangles sorted according to their geodesic distance to a starting vertex.|

Computing exact geodesic distances on 3D triangle meshes is a challenging problem. Traditional shortest path algorithms, such as Dijkstra's approximation algorithm, only traverses along the edges of the mesh, which can cause significant accurancy discrepancies in sparsely sampled meshes. The Mesh Geodesy library is an attempt to provide a solution to this problem by calculating accurate (albeit with floating point approximations) distances and paths for a given vertex to all other vertices.

<img src="pictures/cow_path.png " alt="alt text" width="300" height="250"> <img src="pictures/dancer_path.png " alt="alt text" width="300" height="250"> <img src="pictures/hand_path.png " alt="alt text" width="300" height="250"> <img src="pictures/nicolo_path.png " alt="alt text" width="300" height="250">

This library implements the exact algorithim presented in the paper "Fast Exact and Approximate Geodesics on Meshes" by Hoppe and team. [link][1]

[1]: https://hhoppe.com/geodesics.pdf

## Content

This repository contains several main components:

1. **Mesh Geodesy Library (mblib)**: This is a self-contained library that provides a set of APIs for calculating geodesic distance and path on 3D triangle meshes. It includes APIs for loading a mesh and computing geodesic distances and paths for any given vertex of the mesh.
2. **Mesh Geodesy GUI Application (mgapp)**: This is a GUI-based tool built on FLTK that allows for interactive visualization and manipulation of 3D triangle meshes. The application includes features such as selecting vertices, changing textures, and more.
3. **Mesh Geodesy Command Line (mgcmd)**: This is a command line application that allows for the computation of geodesic distance and path via the command line. The input 3D model and output distance information can all be specified as call arguments.
4. **Mesh Models (models)**: This directory contains testing 3D models that have been cleaned up to demonstrate the properties of a manifold.

## Installation & Requirements

The CMake builder will try to find the following dependency packages/libraries. If not found, error messages with proper installation instructions for that library will be displayed. For best installation experience, the following packages can be installed in advance, to avoid interruptions in the CMake build:
* FLTK - 
```$ sudo apt-get install libfltk1.3-dev```
* OpenGL - 
```$ sudo apt-get install libgl1-mesa-dev```
* Xft -
```$ sudo apt-get install libxft-dev```
* fontconfig - 
```$ sudo apt-get install libfontconfig1-dev```
* Xrender -
```$ sudo apt-get install libxrender-dev```
* Xfixes - 
```$ sudo apt-get install libxfixes-dev```
* Xcursor -
```$ sudo apt-get install libxcursor-dev```
* Xinerama - 
```$ sudo apt-get install libxinerama-dev```
* glog - 
```$ sudo apt-get install libgoogle-glog-dev```
* libpng - 
```$ sudo apt-get install libpng-dev```


Finally, the installation of Mesh Geodesy package can be done by the typical
```
$ mkdir build
$ cd build
$ cmake ..
$ make
```
To launch the GUI mgguiapp, run
```
$ mgapp/mgapp
```
To use the command line function, use
```
$ mgcmd/mgcmd <options>
```
For more details on how to use GUI or command line, please refer to [examples](./examples.md)

### Remote Desktop Connection instructions for Windows users

In order to launch and visualize the GUI windows, if you are under a Windows remote system, the recommendation is to install Remote Desktop Connection by following the instructions below:
```
sudo apt-get update
sudo apt-get install xrdp
```
Once the installation is complete, start the xrdp service by running the following
```
sudo systemctl enable xrdp
sudo systemctl start xrdp
```
Now, open the Remote Desktop Connection application on your Windows computer; Enter the IP address of your system in the "Computer" field and click "Connect"; Follow the instructions to enter username and passward. You should be connected to your Linux system.

### Tested environments

This project was mainly developed and tested under an Ubuntu system: Ubuntu 22.04.2 LTS; with gcc v11.3.0; cmake v3.26.0. Due to development capacity, it has not been extendedly under other environments. Contributions are welcome!

## Known issues

It is know that the implementation requires a manifold 3d mesh. Meaning that the input 3D mesh should satisfy the following criteria:
* Each edge in the mesh is shared by exactly two faces. 
* The mesh should be orientable, meaning that it is possible to consistently assign a normal direction to each face such that the faces form a continuous surface.
* The mesh should be edge-manifold, meaning that the 1-ring neighborhood of each vertex (the set of vertices connected to the vertex by an edge) should form a single connected component, and should not form any "holes" or "handles." This ensures that the local neighborhood of each vertex resembles a flat disk, rather than a more complex topology.


## Acknoledgemeent

This library uses some code from LibQSlim, a library for triangle mesh simplification created by my thesis advisor, Michael Garland. I tried my best to package those files in an independent subdirectory called `gfx/`, but some files such as arcball, GUI were spread in `mgapp/`. I would like to thank Prof. Garland for his invaluable teaching and unintentional contributions.

I would also like to acknowledge the help of ChatGPT for providing advice on using CMake in this project. It made possible to deliver this project (even not yet complete) in such short period of time.


## Personal Notes

Geodesic distance and path calculation was a dependency for my PhD thesis. At the time, there was no available open source that implements the function. Thus, I spent several months to write the code and get all corner cases right. 

While I acknowledge that there may be newer and better solutions available today, I still hope that this library can be useful to others. I still remember the countless hours of debugging and troubleshooting that went into making this library work. If this library can save someone a few weeks of work, or a few sleepless nights, then I consider it a success.
