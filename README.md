# PE-Planner

**PE-Planner** is a performance-enhanced quadrotor motion planner for autonomous flight in complex and dynamic environments. It is proposed to significantly improve the performance of speed, safety, and disturbance rejection capability.

__Authors__: Jiaxin Qiu, Qingchen Liu, Jiahu Qin, Dewang Cheng, Yawei Tian and Qichao Ma

Please give us a :star: if this project helps you! We will open source more.

<p align="center">
  <img src="gif/github_video1.gif" width = "400" height = "225"/>
  <br>
  <img src="gif/github_video2.gif" width = "400" height = "225"/>
  <img src="gif/github_video3.gif" width = "400" height = "225"/>
  <img src="gif/github_video4.gif" width = "400" height = "225"/>
  <img src="gif/github_video5.gif" width = "400" height = "225"/>
  <img src="gif/github_video6.gif" width = "400" height = "225"/>
  <img src="gif/github_video7.gif" width = "400" height = "225"/>
</p>

## Table of Contents

* [Installation](#1-installation)
* [Run Simulations](#2-run-simulations)

## 1. Installation
The project is developed on Ubuntu 20.04 (ROS Noetic). In simulations, it uses Gazebo as the simulator and PX4 as flight control software to achieve relatively realistic simulations. To avoid tedious configuration steps, a Docker image containing the simulation environment and PE-Planner is provided. Running the following commands to setup:

```bash
docker pull kuguao/pe-planner-simulation-env:latest
```
If your computer has an NVIDIA device, run
```bash
docker run -it --ipc=host --volume="/tmp/.X11-unix:/tmp/.X11-unix:rw" --gpus all --env="DISPLAY" --env="NVIDIA_DRIVER_CAPABILITIES=all" --env="QT_X11_NO_MITSHM=1" --name="pe-planner" kuguao/pe-planner-simulation-env:latest
```
Otherwise, run
```bash
docker run -it --ipc=host --volume="/tmp/.X11-unix:/tmp/.X11-unix:rw" --env="DISPLAY" --name="pe-planner" kuguao/pe-planner-simulation-env:latest
```

## 2. Run Simulations
### Simulation of the Nominal Case with Static Obstacles
Open three terminals to run commands individually and sequentially.

Terminal 1:
```bash
docker exec -it pe-planner /bin/bash
cd ~/pe-planner-simulation
./startpx4.sh
```
Terminal 2:
```bash
xhost +
docker exec -it pe-planner /bin/bash
cd ~/pe-planner-simulation/PE-Planner
git checkout static_env_simulation
rviz -d ./rviz/rviz.rviz
```
Terminal 3:
```bash
docker exec -it pe-planner /bin/bash
cd ~/pe-planner-simulation/PE-Planner
./build.sh
./build/planner_px4 ./maps/map3.txt 1 1 6.0 3.0 1.0 25
```
The format of the command to run PE-Planner is 
```
./build/planner_px4 [map_file] [enable MPCC (disabled it to use PID for comparison)] [enable GPIO] [max velocity in Kinodynamic searching] [max acceleration in Kinodynamic searching] [mu] [number of tests]
```

### Simulation of the Nominal Case with Static and Dynamic Obstacles
Terminal 1:
```bash
docker exec -it pe-planner /bin/bash
cd ~/pe-planner-simulation
./startpx4.sh
```
Terminal 2:
```bash
xhost +
docker exec -it pe-planner /bin/bash
cd ~/pe-planner-simulation/PE-Planner
git checkout dynamic_env_simulation
rviz -d ./rviz/rviz.rviz
```
Terminal 3:
```bash
docker exec -it pe-planner /bin/bash
cd ~/pe-planner-simulation/PE-Planner
./build.sh
./build/planner_px4 ./maps/map2.txt 1 1 6.0 3.0 1.0 25 0
```
The format of the command to run PE-Planner is 
```
./build/planner_px4 [map_file] [enable MPCC (disabled it to use PID for comparison)] [enable GPIO] [max velocity in Kinodynamic searching] [max acceleration in Kinodynamic searching] [mu] [number of tests] [enable dynamic planning (disable it when enabling MPCC)]
```
### Simulation of the Case with Disturbance
Terminal 1:
```bash
docker exec -it pe-planner /bin/bash
cd ~/pe-planner-simulation
./startpx4.sh
```
Terminal 2:
```bash
xhost +
docker exec -it pe-planner /bin/bash
cd ~/pe-planner-simulation/PE-Planner
git checkout disturbance_simulation
rviz -d ./rviz/rviz.rviz
```
Terminal 3:
```bash
docker exec -it pe-planner /bin/bash
cd ~/pe-planner-simulation/PE-Planner
./build.sh
./build/planner_px4 ./maps/map2.txt 1 1 6.0 3.0 1.0 25 8.48
```
The format of the command to run PE-Planner is 
```
./build/planner_px4 [map_file] [enable MPCC (disabled it to use PID for comparison)] [enable GPIO] [max velocity in Kinodynamic searching] [max acceleration in Kinodynamic searching] [mu] [number of tests] [value of disturbance force]
```
The value of disturbance force given is only used to set the disturbance in Gazebo simulation and is unknown to the planner.
