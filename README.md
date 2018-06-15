# cuDebye Version 1.5
A fast experimental cuda routine for the evaluation of powder diffraction patterns from atomistic models.

### AT THE MOMENT THIS FILE AND THE BRANCH ARE UNDER CONSTRUCTION DO NOT USE THE PROGRAM UNTIL THIS LINE IS DELETED.

This is an open source software, you can use and redistribute everything, **but you are encouraged to cite our publicaction:**
- will be added soon

In case of commercial use, questions, problems and bugs please write me an [E-Mail](mailto:m.s.rudolph@outlook.com).

**As this is an experimental code the usage is at own risk**, we will not take on responsibility for a system crash or anything else.
_Do not use the program until you have read this documentation file completly._

## GENERAL INFORMATIONS

- the source code will be uploaded soon
- I decided to upload all my programming comments despite of the fact, that they are not written in a proper scientific language and probably contain a lot of spelling mistakes (I have never checked these comments after programming)
- programming language is _C++_ and _Cuda_
- the source code was written and compiled for **_Microsoft Windows (x64)_** using _Microsoft Visual Studio 2015_ and the _NVIDIA Cuda Toolkit 8.0_ --> you will need [_Microsoft Visual C++ 2015 Redistributable (x64)_](https://www.microsoft.com/de-de/download/details.aspx?id=48145)
- at least a [_Cuda Compute Capability_](https://de.wikipedia.org/wiki/CUDA#Unterst%C3%BCtzte_GPUs) of 2.1 or higher is recommended
- the program was tested only on the following graphic processing units (GPU): GTX 970, GTX 980TI and GTX 1080
- the program can use only one GPU and the workload will be around 100%
- **if the same GPU is used for display tasks, your screen will freeze during the calculation**
- thus, I use the _Intel's integrated HD graphics_ for my monitors, but you could also use a second dedicated GPU or wait until the calculation is done

## Input File

As input file for the calculation a so called *.deb-file is necessary. If you check the provided example _GammaAl2O3_10nm.deb_ you will recognize a short header (line 1-4) followed by the different atom/ion types (line 5-7) and their positions (line 8-117943). For simplicity the first 8 lines are given in the following:

1. **0 8 1024**
2. **1.541800 10.00 140.00 0.02**
3. **356 0.0001 0**
4. **3**
5. **67392 O2- 1.0000 2.0000 3.990247 ... 0.025429**
6. **33696 Al3+ 0.8333 2.0000 4.132015 ... 0.019397**
7. **16848 Al3+ 1.0000 2.0000 4.132015 ... 0.019397**
8. **0.0000 3.9115 98.2772**

All values/parameters within a line are separated by space.

### First Line
This line contains probably the most confusing informations and specifies the computing options for the _Cuda Device_ (_Cuda_ capable GPU).
**However this is also the most important line, since wrong parameters can initiate a program crash.**
Usually the provided parameters in _GammaAl2O3_10nm.deb_ are safe for most GPUs and should even work on a [GeForce GT 520](http://www.nvidia.de/object/product-geforce-gt-520-de.html). 
Please do not try GPUs with less power.

**The first number** defines the _Cuda Device_ which is used for the calculation. If you have only one (the onboard GPU, e.g. does not count, as it is not _Cuda_ capable) this number is 0. 
Otherwise each installed device has an integer number 0,1,2... choose the one which is most powerfull.
For checking your PC you can use [_GPU-Z_](https://www.techpowerup.com/gpuz/).

**The second and the third number** specifies the block and the grid size for the calculation.
In the shown example the block size is **8 x 8 = 64 Threads** and the grid size is **1024 x 1024 = 1e6 Blocks = 67e6 Threads** (approximately). 
In other words per calculation step around 67 million distances are calculated, if the number of interatomic distances is higher the program will automatically run enougth steps. 
Increasing the number of distances per step accelerates the whole calculation process, but increases the time necessary for one step.
The latter can be a problem as GPUs are usually reseted by _Microsoft Windows_ after a time of 2 seconds without any respond, to prevent a full display freezing.
Unfortunately this leads to the termination of _cuDebye_.
The faster your GPU, the larger is the number of distances you can calculate within these 2 seconds.

So you should run the program (see **Running cuDebye**) first with the parameters in the example file _GammaAl2O3_10nm_ and check the following displayed **CUDA PROPERTIES**:
- Max Threads per Block --> e.g. 1024
- Threads per Warp --> e.g. 32
- Max Blocks per Grid in Y-Dimension (the one with the lowest value) --> e.g. 65535
Afterwards you can increase the block and grid size carefully under the consideration of the following hints:

**Second Number (Block Size)**
The maximum block size (which can be used) is the square root of maximum threads per block, in our example it is 32.
Furthermore, it is a advantage but not a condition that the numper of threads per block is a multiple of 2 (binary system) and 32 (threads per warp).
It is also an advantage if then numer of threads is not at maximum.

**Third Number (Grid Size)**
The maximum 

First I recommend to run this program 
The maximum size of threads per block and blocks per grids depends on your _Cuda Device_
In other words the GPU calculates 
[_NVIDIA Nsight_](http://www.nvidia.com/object/nsight.html)

