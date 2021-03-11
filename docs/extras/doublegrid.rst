Using the CFM's doublegrid functionality
----------------------------------------

The CFM's numerical scheme is Lagrangian, and each snowfall event gets a new layer (node) added on (and one is removed from the bottom). This creates a lot of nodes (slow!) when time steps are small. So, we developed a 'regrid' or 'doublegrid' function, which combines some specified number of nodes at a specified depth to reduce the number of nodes. Then, to improve this further, the is an additional node combining depth to make an even-lower resolution grid. 

Here is Vincent Verjans' description of this function:

This new scheme is very similar to the original regrid scheme. The difference is that it uses a coarser resolution below grid2, in a grid that I called grid22. So now, you can have a transition in vertical resolution at two different depths. We still have grid1 at high resolution, then grid2 at low resolution and then grid22 at very low resolution. Below grid22, there is the grid23 which is again the same resolution as grid2. As before, the grid3 provides the stock of layers to be removed at each accumulation event. We have thus 5 different grids now:
grid1: high resolution determined by accumulation events
grid2: low resolution
grid22: very low resolution
grid23: low resolution
grid3: high resolution

First, we proceed to an initial gridding in firn_density_spin. This splits the model domain into grid1 - grid2 - grid22 - grid23 - grid3.
The transition depth and the number of nodes to merge from grid1 to grid2 are given by the json entries "grid1bottom" and "nodestocombine" as before. In addition, you can now include in the json the entries "grid2bottom" and "multnodestocombine". The former defines the transition depth from grid2 to grid22 and the latter determines how many nodes of grid2 are merged in a node of grid22. So, 1 node of grid22 is made of ("nodestocombine"x"multnodestocombine") nodes of the grid1. The grid23 nodes are formed of the same number of initial nodes merged together as the grid2 nodes. And the grid3 nodes are formed of a single initial node (as before).
Note that by setting multnodestocombine to 0, the regrid scheme will work exactly as before and there won't be any grid22 and grid23.

Second, in the firn_density_nospin run, regrid22 is called every time the stock of grid3 nodes is empty. At this point, we merge nodes of grid1 into a grid2 node and we split a grid23 node into grid3. But when the stock of grid23 layers is empty, this does not work anymore. So, if this is the case, we merge nodes of grid2 into a node of grid22 and we split a node of grid22 into nodes for grid23. After that, we can again split a grid23 node into nodes for grid3.
There is a good reason we use the grid23 to form grid3 nodes and we don't split a grid22 node into grid3 nodes. The grid22 nodes are very thick. If we split a single grid22 node into "nodestocombine" nodes for grid3, the grid3 nodes would be thicker than the grid1 nodes. On the long term, that would cause a thinning of the total CFM domain because every time a grid1 node is accumulated, a grid3 node is removed. By using grid23, we effectively split a single grid22 node into ("nodestocombine"x"multnodestocombine") nodes for the grid3.