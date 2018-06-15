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
- the source code was written and compiled for **Microsoft Windows (x64)** using _Microsoft Visual Studio 2015_ and the _NVIDIA Cuda Toolkit 8.0_ --> you will need [_Microsoft Visual C++ 2015 Redistributable (x64)_](https://www.microsoft.com/de-de/download/details.aspx?id=48145)
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
Usually the provided parameters in _GammaAl2O3_10nm.deb_ are safe for most GPUs 
The first number defines the _Cuda Device_ which is used for the calculation. If you have only one (the onboard GPU, e.g. does not count, as it is not _Cuda_ capable) this number is 0. 
Otherwise each installed device has an integer number 0,1,2... choose the one which is most powerfull.
For checking your PC you can use [_GPU-Z_](https://www.techpowerup.com/gpuz/).
The second and the third number specifies the block and the grid size for the calculation.
In the example the block size is **16 x 16 = 256 Threads** and the grid size is **4096 x 4096 = 16777216 Blocks = 4.294967296e9 Threads**.
The maximum size of threads per block and blocks per grids depends on your _Cuda Device_
In other words the GPU calculates 
[_NVIDIA Nsight_](http://www.nvidia.com/object/nsight.html)

