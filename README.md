# Project name: happy-tree-friends
## General Guidelines and Requirements
-  Must parallelize your program by at least one of the following methods
    -  OpenMP (minimum requirement)
    -  MPI (get bonus point)
    -  GPU (get extra bonus point)
-  Must provide convincing demonstration of the accuracy and performance of your program
-  Must use GitHub for collaborative development
    -  Tutorial: https://gitbook.tw/
    -  Use the Fork→Pull→Push→Pull Request→Merge workflow
    -  Do NOT just upload the final code → must keep the development history on GitHub
-  Students per group: 2
    -  Inform the TA before April 26
-  Final presentation: June 7
    -  25 mins presentation + 10 mins questions
    -  All students are encouraged to ask ANY questions
-  Bonus points: Surprise Me!!!
## Tree Algorithm
- Implement a tree algorithm to compute gravity in 2D
  - Compare the performance and accuracy with different cell-opening criteria
  - Measure the performance scaling (i.e., wall-clock time vs. number of particles)
- Bonus points
  - You get bonus points automatically since this topic is not covered by our lecture
  - Extend to 3D
  - Orbit integration
  - Periodic boundary condition
- Reference
  - The cosmological simulation code GADGET-2:
    https://wwwmpa.mpa-garching.mpg.de/gadget/gadget2-paper.pdf
  - A hierarchical O(N log N) force-calculation algorithm:
    https://www.nature.com/articles/324446a0


## Meeting records
- https://docs.google.com/document/d/1dTHicoDv07S-JrnXMipeCdmDJ7Y_S8J5LuA3QogfXHE/edit?pli=1

## Other Reference
- Barnes-Hut-Simulator
  - GitHub code: https://github.com/beltoforion/Barnes-Hut-Simulator
  - Article: https://beltoforion.de/en/barnes-hut-galaxy-simulator/

# Installation
```
git clone https://github.com/tim20010803/happy-tree-friends.git
```
<p align="right">(<a href="#readme-top">back to top</a>)</p>

# Implementation

1. 
```
cd initial_condition
```

2. Produce initial data.
If you want to use the data we prepared, skip this step and proceed to the next one.
```
g++ -o produce_data_once produce_data_once.cpp
cd output
./produce_data_once
```
Choose the initial condition you like and copy.

The number indicates the particle number.

"Uni" refers to a uniform particle density.

"Nonuni" refers to an exponential distribution of particle density.For example:
For example:
```
cp one_step_data.csv_Nonuni_1000.csv ../
cd ../../code/
```

3. Use the .sh files to compile and run simulations according to the algorithm you want to use.

"Nontree" refers to the normal algorithm.

"Tree" refers to the algorithm that utilizes a tree structure.

"OMP" stands for Open_MP, indicating the use of parallel processing.

"One step" means performing the calculation without orbit integration.

Before that, modify the particle number in the corresponding main file according to the initial condition you want to use.
For example:
```
vim mainTree_one_step
```
find
```
int pN =1000;//choose particle number 
```
Then, compile and execute.
```
./compileTreeOMP_one_step.sh
cd output
./mainTreeOMP_one_step
```

4.Copy the results and plot.
If the tree algorithm is used, copy the mainTree_data.csv file. Otherwise, copy the mainNontree_data.csv file.
For example:
```
cp mainTree_data.csv ../../plot_and_results
```

