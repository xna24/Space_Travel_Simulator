# Relativistic Space Travel Simulation
## Overview
This Python project simulates relativistic space travel, allowing users to explore the effects of traveling close 
to the speed of light. By configuring various parameters, users can visualize the journey to distant stars, including 
relativistic effects on time, distance and the apparence of stars. 

The basic set-up of the journey is this: we travel along a straightline to the destination star with constant (proper) acceleration
of 9.8 $m\cdot s^{-2}$
(See document.pdf for a complete
introduction of the physics and astronomy behind the project)

## Features
1. Configurable settings via `settings.ini` file
2. Supports Hipparcos; Hipparcos-Tycho or customized star catalogs
3. Adjustable parameters for star plotting (size based on magnitudes)
4. Output in customizable pixel dimensions
5. Option to generate series of images for animation

## Installation
To run this project, you need Python 3 installed on your machine.
1. Clone the repository:
```bash
git clone https://github.com/xna24/Space_Travel_Simulator.git
cd Space_Travel_Simulator
```
2. Installed required packages
```bash
pip install -r requirements.txt
```
3. Configure your settings in the `settings.ini` file

## Usage
You can run the simulation using command line. Here are two main ways to execute the script:
1. Using Default Settings:
```bash
python3 main.py
```
This will use the parameters defined in your `settings.ini` file.
2. Specifying Trip Ratio:
```bash
python3 main.py [trip_ratio]
```
Replace `[trip_ratio]` with a floating-point number between 0 and 1 (not too close to these two numbers), indicating how much of the trip (in proper time) is complete (e.g., 0.5 for halfway).
