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
- I decided to upload all my programming comments, despite of the fact that they probably contain a lot of 
- programming language is _C++_ and _Cuda_
- the source code was written and compiled for **Microsoft Windows (x64)** using _Microsoft Visual Studio 2015_ and the _NVidia Cuda Toolkit 8.0_ --> you will need [_Microsoft Visual C++ 2015 Redistributable (x64)_](https://www.microsoft.com/de-de/download/details.aspx?id=48145)
- at least a [_Cuda Compute Capability_](https://de.wikipedia.org/wiki/CUDA#Unterst%C3%BCtzte_GPUs) of 2.1 or higher is recommended
- the program was tested only on the following graphic processing units (GPU): GTX 970, GTX 980TI and GTX 1080
- the program can use only one GPU and the workload will be around 100%
- **if the same GPU is used for display tasks, your screen will freeze during the calculation**
- thus, I use the _Intel's integrated HD graphics_ for my monitors, but you could also use a second dedicated GPU or wait until the calculation is done

## Input File

As input file for the calculation a so called *.deb-file is necessary. If you check the provided example **GammaAl2O3_10nm.deb** you will recognize a short header (line 1-4) followed by the different atom/ion types (line 5-7) and their positions (line 8-117943)

###First Line

