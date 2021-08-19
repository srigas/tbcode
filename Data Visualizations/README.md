This folder contains Python files used for the visualization of the data acquired from various routines within the programs. Each Python file has a default file to read data from, named pythonfilename_def_data, so that the user can remember what kind of inputs each program requires.

### 📈 2dhist.py

This Python file creates a 2D Histogram where each bin is painted according to the intensity corresponding to that bin. Used in order to visualize the spectral density on and around the site hosting a Majorana Quasiparticle.

### 📈 ebdos.py

This Python file reads the output of the EBDOS.f90 code and plots two Spectral Density plots side by side: the first one has a uniform colormap, while the second one has a non-uniform one, for better resolution.

### 📈 diff_chain_openings.py

This Python file creates 6 blocks of diagrams in a 3x2 grid, where the left part corresponds to EBDOS diagrams, while the right part corresponds to simple line plots.

### 📈 many_atom_ebdos.py

This Python file creates 3 ebdos diagrams, each for a different number of chain lengths. Each diagram has an inset where the topological region has been zoomed in.

### 📈 diff_helix_angles.py

This Python file creates 4 ebdos diagrams, each for a different value of the helix's rotation angle.

### 📈 LvsD.py

This Python file creates two plots side by side, with each plot showing three separate lines and a legend. The example shows how the order parameter varies with the Lambda parameter for a BCC, SC and FCC structure.

### 📈 double_dos.py

This Python file creates two DoS plots side by side, in order to compare the situation where there is no impurity present to the one where a magnetic impurity is present.

### 📈 slab_friedel.py

This Python file creates two line plots, one with an inset with a zoomed region.