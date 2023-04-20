# Project name: happy-tree-friends
## General Guidelines and Requirements
1. Must parallelize your program by at least one of the following methods
  a. OpenMP (minimum requirement)
  b. MPI (get bonus point)
  c. GPU (get extra bonus point)
2. Must provide convincing demonstration of the accuracy and performance of your program
3. Must use GitHub for collaborative development
  a. Tutorial: https://gitbook.tw/
  b. Use the Fork→Pull→Push→Pull Request→Merge workflow
  c. Do NOT just upload the final code → must keep the development history on GitHub
4. Students per group: 2
  a. Inform the TA before April 26
5. Final presentation: June 7
  a. 25 mins presentation + 10 mins questions
  b. All students are encouraged to ask ANY questions
6. Bonus points: Surprise Me!!!
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
