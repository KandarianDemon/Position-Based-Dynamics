### Readme

Attempt to implement Position-based-dynamics as described by Matthias Mueller et. al. 
https://matthias-research.github.io/pages/publications/posBasedDyn.pdf


### Overview

### Mesh Generator
Static class with functions to procedurally generate basic Meshes, like Spheres, Cubes, Octaahedra, Icosahedra and Icospheres. 
Used for the generation of a procedural worm-mesh

### Softbody
Component that initializes Softbody Constraints and runs the simulation

### Constrains

Implementation of different types of constraints. Currently implemented Distance and Bending constraints. Work in progress: Self-collision constraint and Volume-constraint for over-pressurized cloth balloons.

### Test Mesh
Just a test class that combines procedural meshes and softbody simulation.