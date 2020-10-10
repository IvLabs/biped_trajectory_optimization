# Biped Trajectory Optimization
## NOTE : This project is still in developement
- [Biped Trajectory Optimization](#biped-trajectory-optimization)
  * [Dynamic Walking on sinusoidal terrain](#dynamic-walking-on-sinusoidal-terrain)
  * [Dynamic Walking on staired terrain](#dynamic-walking-on-staired-terrain)
  * [Dynamic Walking on sloped terrain](#dynamic-walking-on-sloped-terrain)
  * [Dynamic Walking on flat terrain](#dynamic-walking-on-flat-terrain)
  * [Gait Generation for single step](#gait-generation-for-single-step)
    + [using CasADi library in python](#using-casadi-library-in-python)
  * [Trajectory Optimization on some basic systems](#trajectory-optimization-on-some-basic-systems)
    + [cartpole on python using CasADi](#cartpole-on-python-using-casadi)
    + [simple pendulum](#simple-pendulum)
    + [cartpole on C++](#cartpole-on-c)
  * [Passive Walking of 2-link bipedal system](#passive-walking-of-2-link-bipedal-system)

## Dynamic Walking on sinusoidal terrain
### Human gait
![](five_link_path_generation/uneven_terrain/results/sin_walk_10.gif)
![](five_link_path_generation/uneven_terrain/results/sin_walk_10.png) 
### Ostrich gait
![](five_link_path_generation/uneven_terrain/results/osin_walk_10.gif)
![](five_link_path_generation/uneven_terrain/results/osin_walk_10.png) 
## Dynamic Walking on staired terrain
### Human gait
![](five_link_path_generation/uneven_terrain/results/stairs_walk_10.gif)
![](five_link_path_generation/uneven_terrain/results/stairs_walk_10.png) 
![](five_link_path_generation/uneven_terrain/results/stairs_down_walk_10.gif)
### Ostrich gait
![](five_link_path_generation/uneven_terrain/results/ostairs_walk_10.gif)
![](five_link_path_generation/uneven_terrain/results/ostairs_walk_10.png) 
![](five_link_path_generation/uneven_terrain/results/ostairs_down_walk_10.gif)
## Dynamic Walking on sloped terrain
### Human gait
![](five_link_path_generation/uneven_terrain/results/slope_walk_10.gif)
![](five_link_path_generation/uneven_terrain/results/slope_walk_10.png) 
### Ostrich gait
![](five_link_path_generation/uneven_terrain/results/oslope_walk_10.gif)
![](five_link_path_generation/uneven_terrain/results/oslope_walk_10.png) 
## Dynamic Walking on flat terrain
### Human gait
![](five_link_path_generation/uneven_terrain/results/flat_walk_10.gif)
![](five_link_path_generation/uneven_terrain/results/flat_walk_10.png) 
### Ostrich gait
![](five_link_path_generation/uneven_terrain/results/oflat_walk_10.gif)
![](five_link_path_generation/uneven_terrain/results/oflat_walk_10.png) 

## Gait Generation for single step
### using CasADi library in python
![](five_link_gait_generation/animation2.gif) ![](five_link_gait_generation/graph.png)

## Trajectory Optimization on some basic systems
### cartpole on python using CasADi

![](basic_tasks/catpole_python/cartpole.gif) ![](basic_tasks/catpole_python/Graph.png)

### [simple pendulum](basic_tasks/simple_pendulum.m)

### [cartpole on C++](basic_tasks/cartpole_cpp)

## [Passive Walking of 2-link bipedal system](passive_walker)

