# Examples

## 1. mgcmd

### Step 1 - Find 1:N geodesic on a model

This  example demonstrates how to use the command prompt to perform geodesic calculations using `mgcmd`. To perform the geodesic calculation from vertex 1014 to the rest of the vertices in the input mesh (cow.obj), simply run the following command (from the `build` directory):

````
cd <mesh_geodesy_root>/build/
mgcmd/mgcmd --input_file ../models/cow.obj --source_vertex 1014 --output_file ../models/cow1014.mg 
````
This command will take the input mesh from the `cow.obj` file, compute the geodesic distances from vertex 1014 to all other vertices, and save the results in the cow1014.mg file. You can easily adapt this example to work with different input meshes and source vertices by modifying the command-line arguments accordingly.

### Step2 - Visualize results using `mgapp`

This session will guide you through visualizing the geodesic calculation results using the provided GUI application. First, run the following command to open the GUI window:
```
mgapp/mgapp
```
Once the GUI window is open, you can load the results from the previous step for visualization. To do this, click on "File" menu item then select "Load Mesh Geodesy" option.

<img src="pictures/load_geodesy.png " alt="alt text" width="400" height="180">

After selecting this option, wait for the results to load, and you will see the output displayed in the GUI, as illustrated in the figure below.

<img src="pictures/cow_geodesy.png " alt="alt text" width="500" height="350">

Enable the rendering of the geodesic path by selecting "Draw" menu item, and enable the option "Geodesic Path".

You can interact with the GUI using mouse click and drag controls:
* Mouse drag with button 1: Rotate the object.
* Mouse drag with button 2: Pan the view across the screen.
* Mouse drag with button 3: Zoom in and zoom out.

