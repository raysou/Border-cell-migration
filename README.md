# Simulating border cell migration: An active bead-spring model
This repository contains a FORTRAN code and a sample data file containing the positional information of all the cells over time (for the WT condition). The FORTRAN code first creates an initial organization of the Drosophila egg chamber containing three different cell types: Nurse, Border, and Polar cells. Then it simulates the time evolution of the whole tissue system by numerically solving the equations S6-S9 in the manuscript: 
"**Modeling Collective Cell Migration in Heterogeneous Environments Using Drosophila Border Cells**", Soumyadipta Ray (1), Tuhin Roy (2), Sayan Acharjee (1), Gaurab Ghosh (1), Mohit Prasad (1), and Dipjyoti Das (1)(++), Cell Reports Physical Science, 2026.

(1) Department of Biological Sciences, Indian Institute of Science Education and Research Kolkata, Mohanpur, Nadia, West Bengal, 741246, India.
(2) Interdisciplinary Biophysics Graduate Program, The Ohio State University, Columbus, USA
(++) Corresponding Authors.

Email: raysoumyadipta@gmail.com, dipjyoti.das@iiserkol.ac.in

**Instruction to get the code**: Click on the 'Code' menu (green colour) and from the drop-down menu, choose the option 'Download Zip'. Then extract it to any directory on your desktop.

**Instructions to run the code**: Open the code via any editor, such as 'Document Viewer' or 'Visual Studio Code'. Install the compiler 'gfortran' using the command line: 'sudo install gfortran' (in Linux or WSL in Windows). Then, compile the code as: gfortran -O3 code_name. 
After that, in the same directory, run the code as: ./a.out
All the parameter values are already given inside the code and have been annotated. One can also change those parameter values as well (from Table S1 in the paper). Finally, the code will generate a data file (like the example data file, given here), containing the time evolution of all the cell coordinates. One can make a movie out of it and do other analyses as well. 
**Data structure**: The output data file is structured as follows: Column1: Iteration, Column2: Cell index, Column3: Bead index, Column4: Bead x-coordinate, Column5: Bead y-coordinate, Column6: Migration force x-component value, Column6: Migration force y-component value. Note that cell indices 1 to 6 represent the six nurse cells, 7 to 12 represent the six border cells, and 13 to 14 represent the two polar cells.
