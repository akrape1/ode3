# ode3

physx630 odelib
---

To build the ODE library and example programs, simply type `make` in this top level ode3 directory.

Description of example programs:<br>

**RKnTest**: Solves a single 1st order ODE using the single equation RK4 solver and the ODE array solver

**RKnStep**: A basic example of the ODE array solver is applied to projectile motion with a simple model of air resistance, force of air resistance = -kv^2<br>. At each step in the elapsed time and x,y positions are printed.<br>
Optional parameters [default values]<br>
* -v initial_velocity [100] m/s
* -t angle_thera [45] degrees
* -m mass_of_projectile [10] kg
* -k coefficient_of_air_resistance [0.1] kg/m


**RKnDemo**: Solves for projectile motion with a simple model of air resistance, force of air resistance = -kv^2<br>
This program includes graphical output.  Detailed output is saved in TGraph objects in RKnDemo.root.  The file **RKnPlotDemo.py** shows how to access date in the TGraphs and can be used to generate additional plots.<br>
Optional parameters [default values]<br>
* -v initial_velocity [100] m/s
* -t angle_thera [45] degrees
* -m mass_of_projectile [10] kg
* -k coefficient_of_air_resistance [0.1] kg/m

**baseball1**:  Starter template for first baseball problem

**baseball2**:  Starter template for second baseball problem

**baseball_drag.ipynb**: this notebook describes the drag force equations used in the text.

gsl starter code
---

The starter code here (projGSL.cpp) demonstrates very basic usage of the gsl for solving a problem of coupled differential equations. An 8th order R-K solver with fixed step size is used. You are encouraged to try other solvers as you explore the problem. See here for the gsl docs: https://www.gnu.org/software/gsl/doc/html/ode-initval.html

This example solves the 2D projectile motion problem with a simple model for air resistance. After each step, data are stored in ROOT TGraphs, which are then displayed at the conclusion of the calculation.

The gsl provides a number of ODE solvers and a variety of interfaces.  Some of the solvers (not R-K methods) use the Jacobian matrix, which gives the devivative of the function wrt the dependent parameters.  See the gsl examples for details.

Python starter code
---

Two examples are given for using ODE solvers from the scipy.integrate sub-package in Python. In these examples graphs are made using matplotlib.

    Solution (projScPY2.py[ipynb]) using a more modern interface scipy.integrate.solve_ivp. See also: https://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html and https://www.programcreek.com/python/example/119375/scipy.integrate.solve_ivp

    Solution (projScPY.py[ipynb]) using an older interface scipy.integrate.odeint¶ (see comments here: https://docs.scipy.org/doc/scipy/reference/integrate.html).  I do not recommend using this interface any longer.

The notebook versions contain additional comments on using the integrators.

Part 1:

so the instructions weren't very clear to me. I wrote vterm.cpp to do the first part without air resistance and vterm2.cpp do the air resistance, terminal velocity, and v_t vs mass plot. These two scripts are executable by running make and then running the executable like all the other cpp scripts provided. However, these scripts are based off of RKnDemo.cpp which writes its outputs to a ROOT file (vterm.root or vterm2.root respectively). I work in VSCode with the JSROOT extension to look at ROOT files. A lot of my plots however looked terrible or didn't even show up with this extension. So there are two additional scripts: test.cpp and test2.cpp which are used to generate the E vs t plot (test.cpp) and the v_t vs m plot (test2.cpp). These two scripts can be ran directly in ROOT with .x test.cpp or .x test2.cpp. I much prefer running scripts like this compared to the Makefile stuff. The vterm.pdf has some plots and thoughts, made in overleaf.

here's the terminal velocity. RKnDemo.cpp already calculates it with default parameters so I didn't do anything new. 
Final velocity = 27.6553

Part 2:

So the baseball speed is quite high, but the drag may be pretty high which causes that. The math looks right even if the physics seems not great. 48 m/s is 107 mph. A really quick google shows pitches sub 100 mph, so I think the problem is that the model is too harsh when it comes to drag? that's the only thing I can think of that went wrong since I just followed the layout set in the canvas instructions
********************************
(xend,z0,theta0) = (18.500000,1.400000,1.000000)
v_pitch = 48.255908 m/s
********************************

Also I put baseball1og.cpp in the src folder. This is the original baseball1.cpp script before I started messing with everything. Similarly baseball2og.cpp for the last part does the same thing but for baseball2.cpp. My scripts diverged enough from the starters that I felt the need to store the originals. They cannot be run like this since I didn't modify the Makefile to accomodate them, but the original version of this assignment can do that. 

Part 3:

While my code doesn't follow the starter structure again, it still uses the format of running like ./baseball2 -p i where i=0,1,2,3 for the different pitches. Each pitch is plotted and saved as a png. I threw all the plots into one pdf called pitches.pdf as stated in the directions. I made this in overleaf. Also I forgot to copy and paste the coordinates and speeds, but the script does those, so running should give them. woops!

After I finished eveyrthing, I made another folder called "outputs" to store the ROOT files and the graphs. I just updated test and test2 since they call on the ROOT files. Hopefully that doesn't break anything.
