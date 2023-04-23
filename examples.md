# Examples

## 1. Using `mgcmd` to get the distances and paths for a cow

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

Navigate to `<mesh_geodesy_root>/models` directory, and load `cow.obj`. Wait for the 3D model to load and render. Enable the rendering of the geodesic path by selecting "Draw" menu item, and the option "Geodesic Path". Alternatively, to enable the rendering of the distance information, Select "Draw" menu, deselect all drawing options, except "Geodesic Distance". You will see the output displayed in the GUI, as illustrated in the figure below, respectively.

<img src="pictures/cow_geodesy.png " alt="alt text" width="500" height="350">
<img src="pictures/cow_distance.png " alt="alt text" width="500" height="350">

You can interact with the GUI using mouse click and drag controls:
* Mouse drag with button 1: Rotate the object.
* Mouse drag with button 2: Pan the view across the screen.
* Mouse drag with button 3: Zoom in and zoom out.

Use the following keys to control the frequency of the isolines that represents the distance information:
* Up to make the frequency higher, thus denser iso-lines
* Down to make the frequency lower, thus sparser iso-lines.

## 2. Using the `mgapp` to get the distances and paths for `Lucy`
### Step 1 - Load the mesh

This session will guide you through calculating and visualizing the geodesic calculation results using the provided GUI application. First, run the following command to open the GUI window:
```
cd <mesh_geodesy_root>/build/
mgapp/mgapp
```
Click on "File" menu item then select "Load Mesh" option. Navigate to `<mesh_geodesy_root>/models` directory, and load `lucy.obj` file. Wait for the 3D model to load and render. 

### Step 2 - Select the source vertex

In this step, you will select a vertex using mouse click to be the source of the geodesic propagation. To select a vertex, use the following mouse interactions:
* Mouse drag with button 1: Rotate the object.
* Mouse drag with button 2: Pan the view across the screen.
* Mouse drag with button 3: Zoom in and zoom out.
* Shift + mouse click on a mesh vertex to select the vertex. If a vertex is properly selected, it will be rendered with a red sphere around the vertex.

** caveat: Vertex picking in this application can be somewhat challenging due to the requirement for the mouse click position to be extremely close to the actual vertex of the mesh, without any ambiguity caused by surrounding or overlapping vertices. This is particularly relevant for models with densely packed vertices, such as the Lucy model. To successfully pick a vertex, you may need to zoom in significantly. When a vertex is successfully selected, it will be highlighted with a small red sphere, as shown in the image below.

<img src="pictures/lucy_selection.png " alt="alt text" width="500" height="350">

### Step 3 - Computing distances and paths
When the vertex selection is successful, you can click on "Tools" menu item, and then select "Compute Dist and Paths" option. Alternatively, press `Ctrl + D`to start the core geodesic computation. The computation may take a while (a few mins - more details about performance to follow). Wait until it finishes.

### Step 4 - Visualizing results

Use the following keys to control the frequency of the isolines that represents the distance information:
* Up to make the frequency higher, thus denser iso-lines
* Down to make the frequency lower, thus sparser iso-lines.

<img src="pictures/lucy_distance.png " alt="alt text" width="500" height="350">