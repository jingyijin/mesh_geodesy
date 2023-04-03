# Mesh Geodesy


Geodesic distance calculation is a powerful technique in mesh processing that allows measurement of the shortest path between two points on a 3D surface. Unlike Euclidean distance, which assumes a flat space, geodesic distance takes into account the curvature of the surface and provides a metric that is constrained to the geometry of the mesh. This makes it particularly useful for a variety of applications, in computer graphics, computational geometry, physical simulations, and other areas.

![alt text](pictures/lucy3.gif "Progressive rendering of a 3D mesh with triangles sorted according to their geodesic distance to a starting vertex.")

Computing exact geodesic distances on 3D triangle meshes is a challenging problem. Traditional shortest path algorithms, such as Dijkstra's approximation algorithm, only traverses along the edges of the mesh, which can cause significant accurancy discrepancies in sparsely sampled meshes. The Mesh Geodesy library is an attempt to provide a solution to this problem by calculating accurate (albeit with floating point approximations) distances and paths for a given vertex to all other vertices.

<img src="pictures/cow_path.png " alt="alt text" width="300" height="250"> <img src="pictures/dancer_path.png " alt="alt text" width="300" height="250"> <img src="pictures/hand_path.png " alt="alt text" width="300" height="250">

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
* glog - ```sudo apt-get install -y libgoogle-glog-dev```
* libpng - ```sudo apt-get install libpng-dev```

## Acknoledgemeent

This library uses some code from LibQSlim, a library for triangle mesh simplification created by my thesis advisor, Michael Garland. I would like to thank him for his invaluable teaching and unintentional contributions.

I would also like to acknowledge the help of ChatGPT for providing advice on using CMake in this project. It made possible to deliver this project (even not yet complete) in such short period of time.

[2]: https://github.com/fltk/fltk/blob/master/README.txt

## Personal Notes
This library was developed as part of my PhD thesis. It took me several months to implement it handling different corner cases, and I still remember the countless hours of debugging and troubleshooting that went into making this library work.

While I acknowledge that there may be newer and better solutions available today, I still hope that this library can be useful to others. My motivation for contributing to the open source community is to save researchers and developers  time and effort. If this library can save someone a few weeks of work, or a few sleepless nights, then I consider it a success.
