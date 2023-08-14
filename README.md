
# RoboGator

RoboGator is a robotic navigation project that uses Player/Stage for simulation. This repository contains configuration files, world descriptions, and a simple client for controlling the robot in a simulated environment.

## Prerequisites

- Player/Stage installed on your machine.
- A C++ compiler (like g++).

## Repository Structure

- `CMakeLists.txt`: The CMake configuration file for building the project.
- `formatted_binary_floor.bits`: A file related to the world description.
- `map.inc`: Configuration for the map in the simulation.
- `player_simple.cfg`: Player configuration file.
- `robot.inc`: Robot configuration for the simulation.
- `sick.inc`: Configuration for the SICK laser range-finder sensor.
- `simple.world`: World description for the Player/Stage simulation.
- `simple_client.cpp`: The main client code for controlling the robot.

## Setup & Running

1. **Clone the repository:**
    ```
    git clone https://github.com/EliyaHouri/RoboGator.git
    ```
2. **Navigate to the project directory:**
    ```
    cd RoboGator
    ```
3. **Build the project:**
    ```
    mkdir build
    cd build
    cmake ..
    make
    ```
4. **Run the simulation:**
    In one terminal, in the project dir, launch the Player with the configuration file:
    ```
    player player_simple.cfg
    ```
    In another terminal, in the build dir, run the simple client to control the robot:
    ```
    ./simple_client
    ```

## Contributions

Feel free to fork the project, make some updates, and submit pull requests. Contributions are welcome!
