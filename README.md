# Mesh Geodesy

Geodesic distance calculation is a powerful technique in mesh processing that allows measurement of the shortest path between two points on a 3D surface. Unlike Euclidean distance, which assumes a flat space, geodesic distance takes into account the curvature of the surface and provides a metric that is constrained to the geometry of the mesh. This makes it particularly useful for a variety of applications, such as computer graphics, computational geometry, and physical simulations.

Computing exact geodesic distances on 3D triangle meshes is a challenging problem. Traditional shortest path algorithms, such as Dijkstra's algorithm,  only traverses along the edges of the mesh, which can cause significant discrepancies in sparsely sampled meshes. The Mesh Geodesy library is an attempt to provide a solution to this problem by calculating accurate (albeit with floating point approximations) distances and paths for a given vertex to all other vertices.

This library implements the exact algorithim presented in the paper "Fast Exact and Approximate Geodesics on Meshes" by Hoppe and team. [link][1]

[1]: https://hhoppe.com/geodesics.pdf

## Content

This repository consists of two main components:

### Mesh Geodesy library
mblib: It is a self contained library that provides a set of APIs for calculating geodesic distance and path on 3D triangle meshes. It comes with APIs that allows loading of a mesh and compute geodesic distances and paths for any given vertex of the mesh. 

### Mesh Geodesy Application
mgapp: It is a GUI-based tool built on FLTK that allows interactive visualization and manipulation of 3D triangle meshes. The application provides a few more features, such as selecting vertices, changing textures, and more.

## Installation Requirements

This project was tested under a Ubuntu system. It requires installation of the following dependencies
* FLTK - Installation instructions available through [link][2]
[2]: https://github.com/fltk/fltk/blob/master/README.txt
* glog - ```sudo apt-get install -y libgoogle-glog-dev```
* libpng - ```sudo apt-get install libpng-dev```

## Acknoledgemeent

This library uses some code from LibQSlim, a library for triangle mesh simplification created by my thesis advisor, Michael Garland. I would like to thank him for his invaluable teaching and unintentional contributions.

I would also like to acknowledge the help of ChatGPT for providing advice on using CMake in this project. It made possible to deliver this project (even not yet complete) in such short period of time.
