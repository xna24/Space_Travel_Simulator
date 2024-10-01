# Relativistic Space Travel Simulation
## Overview
This Python project simulates relativistic space travel, allowing users to explore the effects of traveling close 
to the speed of light. By configuring various parameters, users can visualize the journey to distant stars, including 
relativistic effects on time, distance and the apparence of stars. 

The basic set-up of the journey is this: we travel along a straightline to the destination star with constant (proper) acceleration
of $g=$ 9.8 $m\cdot s^{-2}$ for the first half and $-g$ for the second half to slow down before arrival. See `document.pdf` 
for a brief introduction and computations involved. A general rule of thumb is for a trip to a star $L$ lightyears away, the
maximum Lorentz factor is $L/(2\text{ly})$. The proper time experienced by the spaceship crew is in general slower
than time in rest frame. Trip ratio is defined to be a number between 0 and 1 indicating how much
the trip is completed in terms of proper time.

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
alternatively
```bash
pip install numpy matplotlib pandas
```
3. (Optional) Uncompress the Tycho catalog file `HIP_TYC.csv.gz` if using `HIP_TYC` in `StarCatalog` option in `settings.ini`. Be sure to put `HIP_TYC.csv` in the same folder as `main.py`.
4. Configure your settings in the `settings.ini` file

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

## Configuration
The `settings.ini` file allows you to customize various parameters for the simulation. Below is an explanation of the key sections and their options:
### Catalog Section
* StarCatalog: specifies which catalog to use. Options include `HIP_TYC` for the Hipparcos-Tycho catalog (combining stars in both catalogs and uses the more precise parameters in Hipparcos when possible, a total of more than $10^6$ stars); `HIP` for just the Hipparcos catalog (118,218 stars); or `custom` for a user-defined catalog.
* StarCatalogFileName: path to the custom star catalog file, if using a custom catalog
* StarCatalogRAColumn, StarCatalogDEColumn, StarCatalogParallaxColumn, StarCatalogVmagColumn, StarCatalogBVColumn: these parameters define the name of the columns in the star catalog that contain the right ascension (RA), declination (DE), parallax, visual magnitude (Vmag), and B-V color index respectively.

### Destination Section
* HIP: the index in the Hipparcos catalog for the destination star. For custom destinations, set this to `0`.
* Name: the name of the destination star. If left blank and not a custom destination, the HIP number will be printed instead
* RAdeg: right ascension in degree of the destination star (only for custom destination)
* DEdeg: declination in degree of the destination star (only for custom destination)
* parallax: parallax in miliarcseconds (mas) of the destination star (only for custom destination)
To convert between parallax and distance: D=1000/parallax where parallax is in mas and distance D is in parsec. 1 parsec = 3.62 ly

### Sky Map Section
* Width, Height: dimension of the output image in pixels. (These are only a guideline as matplotlib adds margins which I did not find a way to control. The default value left in `settings.ini` results in a roughly 1920x1080 image)
* LaTeX: set to `yes` or `no` to determin whether to render text using LaTeX for better visual quality.
* SizeMultiplier: a multiplier for adjusting the size of stars based on their magnitudes.
* MaxSize: maximum size for the stars in the output image
* SizePowerLaw: star size is determined by this function (const) $\cdot b^{m_V}$ with $m_V$ the V magnitude, (const) is `SizeMultiplier` above and this value is $\log_{10}b$. The default value is 2.5, so that the area of the star image is proportional to the actual intensity but it might result in too few stars visible.
* AzimuthSpan: the azimuthal span of the view, in degrees
* AltitudeSpan: the altitude span of the view, in degrees
* VmagCutoff: the visual magnitude cutoff, stars fainter than this magnitude will not be displayed

### Animation Section
* Animated: set to `yes` or `no` to determine if a series of images will be produced for animation. If not animated, a single file `sky_map.png` will be generated.
* TripStar, TripEnd, TripStep: these will produce an equally spaced trip-ratios to plot
* FrameFormat: the file extension of frame images
* FrameFileHead, FrameFileTail: the strings to attach before / after frame number, e.g. if we plot 150 frames and set
```ini
[Animation]
...
FrameFormat = png
FrameFileHead = frame
FrameFileTail = 20241001
```
The resulting files will be `frame_001_20241001.png`, `frame_002_20241001.png` ...

## Example

A trip to HIP 38594. Initialize the parameters (parts omitted ... left as default) 
```ini
[Catalog]
StarCatalog = HIP_TYC
...
[Destination]
HIP = 38594
...
[Sky Map]
Width = 2405
Height = 1377
LaTeX = no
SizeMultiplier = 2
MaxSize = 25
SizePowerLaw = 2.5
AzimuthSpan = 30
AltitudeSpan = 17.1767
VmagCutoff = 16
...
```

In command-line:
```bash
> python3 main.py 0.15
Star catalog loaded.
Tot. star count: 1040740
Destination: HIP 38594 HIP 38594
Distance: 63.752 ly
Right Ascension(deg): 118.546
Declination(deg): -25.304
██████████████████████████████████████████████████| 0.06%, 2.39s
██████████████████████████████████████████████████| 0.12%, 4.87s
██████████████████████████████████████████████████| 0.18%, 7.16s
██████████████████████████████████████████████████| 0.24%, 8.82s
██████████████████████████████████████████████████| 0.30%, 13.94s
██████████████████████████████████████████████████| 0.36%, 19.78s
██████████████████████████████████████████████████| 0.42%, 23.09s
██████████████████████████████████████████████████| 0.49%, 24.78s
```

Here is the output `sky_map.png`:
![`sky_map.png`](https://github.com/xna24/Space_Travel_Simulator/blob/main/sky_map.png?raw=true)
