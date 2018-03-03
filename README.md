# vinaSH
VinaSH: Improved Autodock Vina that accounts for Sigma Hole interactions.

Citation: If you you use vinaSH. Please cite the following manuscript.
"S---O, S---N sulfur bonding interactions in Protein-Ligand Complexes: Empirical Considerations and Scoring Function” Journal of chemical Information and modeling, 2016, 56 (12), pp 2298–2309


VinaSH Installation Instructions:

1)	Download VinaSH from https://github.com/ssirimulla/vinaSH or from http://www.sirimullaresearchgroup.com/software.html
2)	Save the VinaSH directory in file system at a location of your choosing.
Windows
Compatibility
VinaSH is expected to work on Windows XP and newer systems.
Installing
Download the vinaSH.exe from the Windows folder into the desired location
Running
Open the Command Prompt and go to the folder where you downloaded the executable and run vinaSH.exe --help
Linux
Compatibility
VinaSH is expected to work on x86 and compatible 64-bit Linux systems.
Installing
Download the executable file in the linux folder to the desired location
Running
In a terminal window, go to where you downloaded the executable and run 
./vinaSH --help
If you get a permissions denied message, use the command: chmod +x vinaSH
Alternatively, before downloading the executable, you can go to the directory you wish to work and use the command git clone https://github.com/ssirimulla/vinaSH.git
Then navigate to the linux folder and run ./vinaSH --help
Mac
Compatibility
VinaSH is expected to work on Mac OS X 10.6 (Snow Leopard) through 10.10 (Yosemite) Intel machines.
Installing
	Download the executable file in the mac folder to the desired location 
Running
In a terminal window, go to where you downloaded the executable and run 
./vinaSH --help
If you get a permissions denied message, use the command: chmod +x vinaSH
Alternatively, before downloading the executable, you can go to the directory you wish to work and use the command git clone https://github.com/ssirimulla/vinaSH.git
Then navigate to the mac folder and run ./vinaSH --help
Building from source
Step 1: Install a C++ compiler suite
On Windows, you may want to install Visual Studio; on OS X, Xcode; and on Linux, the GCC compiler suite.


Step 2: Install Boost
Install Boost. (Version 1.57.0 was used to compile the official binaries. With other versions, your luck may vary) 

Step 3: Build Vina
If you are using Visual Studio, you may want to look at this tutorial for compiling: https://sites.google.com/site/mkoohim/stories/how-can-compile-autodock-vina-in-visual-studio. For optimal performance, remember to compile using the Release mode.
On OS X and Linux, you may want to navigate to the appropriate build subdirectory, customize the Makefile by setting the paths and the Boost version, and then type: make depend
And then type: make
3)	VinaSH needs supporting visualization programs for preparation and viewing of molecules. These can be downloaded from any compatible 3rd party:
a.	Some examples are given below
MGL Tools (http://mgltools.scripps.edu/downloads)
Pymol (https://www.pymol.org/)
4)	Once VinaSH and associated programs are downloaded you are ready to start preparing files to run
5)	This program is very similar to the base program AutoDock Vina, for any further information needed please refer to the Scripps Research Institute manual for AutoDock Vina.
http://vina.scripps.edu/manual.html
