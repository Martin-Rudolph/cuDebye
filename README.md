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
- I decided to upload all my programming comments despite of the fact, that they are not written in a proper scientific language and contain a lot of spelling mistakes (in fact the english is horrible, but I have no muse to check the comments)
- programming language is _C++_ and _Cuda_
- the source code was written and compiled for **_Microsoft Windows (x64)_** using _Microsoft Visual Studio 2015_ and the _NVIDIA Cuda Toolkit 8.0_ --> you will need [_Microsoft Visual C++ 2015 Redistributable (x64)_](https://www.microsoft.com/de-de/download/details.aspx?id=48145)
- copying the source code can introduce errors, e.g. <<< LaunchKernel >>> to << < LaunchKernel >> >, commonly the compiler can handle that
- at least a [_Cuda Compute Capability_](https://de.wikipedia.org/wiki/CUDA#Unterst%C3%BCtzte_GPUs) of 2.1 or higher is recommended
- the program was tested only on the following graphic processing units (GPU): GTX 970, GTX 980TI and GTX 1080
- the program can use only one GPU and the workload will be around 100%
- **if the same GPU is used for display tasks, your screen will freeze during the calculation**
- thus, I use the _Intel's integrated HD graphics_ for my monitors, but you could also use a second dedicated GPU or wait until the calculation is done


## Input File
As input file for the calculation a so called *.deb-file is necessary. When you check the provided example [_GammaAl2O3_10nm.deb_](Data/GammaAl2O3_10nm.deb) you will recognize a short header (line 1-4) followed by the different atom/ion types (line 5-7) and their positions (line 8-117943). For simplicity the first 8 lines are given in the following:

![Input File](/images/header.PNG)

All values/parameters within a line are separated by space.

### First Line
This line contains probably the most confusing informations and specifies the computing options for the _Cuda Device_ (_Cuda_ capable GPU).
**Unfortunately, this is also the most important line, since wrong parameters can initiate a program crash.**
Usually the provided parameters in [_GammaAl2O3_10nm.deb_](Data/GammaAl2O3_10nm.deb) are safe for most GPUs and should even work on a [GeForce GT 520](http://www.nvidia.de/object/product-geforce-gt-520-de.html). 
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

So you should run the program (see **Running cuDebye**) first with the parameters in the example file [_GammaAl2O3_10nm.deb_](Data/GammaAl2O3_10nm.deb) and check the following displayed **CUDA PROPERTIES**:

![Console](/images/cmd_CudaProp.PNG)

Afterwards you can increase the block and grid size carefully under the consideration of the following hints:

**Choosing the block size:**
The maximum block size (which can be used) is the square root of maximum threads per block, in our example it is 32.
Furthermore, it is an advantage but not a condition, that the number of threads per block is a multiple of 2 (binary system) and 32 (threads per warp).
It is also an advantage if then number of threads is not at maximum.

**Choosing the grid size:**
The maximum possible grid size is exactly the smallest grid dimension, in our example it is 65636 as given for the y-dimension. Here it is also beneficial when the number is a multiple of 2.
Usually for grid sizes exceeding 8192 the speed improvement is negligible and the monitioring of the calculation progress is delayed.

**Modify the GPU time out:**
As already mentioned the GPU time out is about 2 seconds. Installing [_NVIDIA Nsight_](http://www.nvidia.com/object/nsight.html) allows to change the time out for larger block and grid sizes. You can also disable the time out, but I would not recommend this.

![Nsight 5.2](/images/Nsight1.PNG)

![Nsight Monitor Options](/images/Nsight2.PNG)

**I prefer as block size 16 and as grid size 4096 or 8192 (GTX 980TI), the time out was set to 60 seconds.**

### Second Line
This line contains the information about the diffractogram, which is calculated.

**The first number** is the wavelength, in our case it is 1.5418 Angstrom (CuKa).

**The following numbers** are the minimum diffraction angle, the maximum diffraction angle and the step width, respectively.
All values are given in degree.
The lower and upper limits are 0 and 180 degree.
However, it is recommended to use typical diffractometer values like 10-140 degree with a step size of 0.02 degree.

### Third Line
This line contains the informations for the distance calculation and also specifies which output files are generated.

**The first and the second number** are the maximum occuring interatomic distance and the step width, respectively. These values are used to allocate an integer array for counting the distances. In our case this array comprises **356/0.0001 = 3.6e6 numbers = 30 MB** (approximately). The memory location is associated with a discrete distance, when a distance is calculated the counter for this distance is increased by one. **Choosing a maximum distance lower than the actual maximum distance between the atom positions will result probably in a program crash**, because non reserved memory locations are overwritten. Therefore, it is better to choose always a distance larger than the actual cluster size (around 5%). Increasing the step width reduces also the memory necessary but reduces also the calculation accuracy significantly. Decreasing the step width do not improve the pattern quality significantly, as only single precision operations are used for the distance calculation. The value 0.0001 is recommended. **If you have large clusters keep an eye on your GPU memory**. Start with small clusters and check the memory consumption. Commonly with 3GB or more dedicated GPU memory there should be no problems.

![Free Memory](/images/cmd_Memory.PNG)

**The third number** defines the generated output files:
- 0: only a intensity file (debI) is generated, the first column contains the diffraction angles and the second column the calculated intensities
- 1: only distance histogram files are generated for each subgroup combination, the first column contains the distances in Angstrom and the second column the frequencies (e.g. the file for the interatomic distances between the subgroups m and n is terminated by debxmxn)
- 2: intensity file and distance histograms are created

Writing distance histograms will slow down the routine.
In our example we have 3 different subgroups, thus **(3²+3)/2 = 6** distance histogram files are generated.

![Output Files](/images/output.PNG)

### Fourth Line
This line contains only the number p of different subgroups.
A subgroup comprises all atoms/ions with the same properties, like charge, occupancy, isotropic temperature factor and atomic scattering factors.

### Fifth Line up to Line 4+p
Every line contains the following information for one subgroup.
The subgroups are given in ascending order.
_cuDeby_ will number the first subgroup with index 0.

**The first number** is the number of atoms/ions in the subgroup.

**The second argument** is the atom/ion name.

**The third and fourth number** are the occupacy and the isotropic temperature factor, respectively.

**The last eleven numbers** are the parameters a1-a5, b1-b5 and c for the atomic scattering factor. A list copied from [_RIETNAN_](http://www.ccp14.ac.uk/ccp/web-mirrors/rietan/fujioizumi/rietan/angle_dispersive/angle_dispersive.html) can be found [here](/Matlab/Tables/atomic_scattering_factors_RIETAN.txt).

![Atomic Scattering Factors](/images/asf.PNG)

For neutron diffraction the scattering factors are constant, set all a and b parameters to zero and use a proper c value.

### Line 5+p and Following
All further lines contain the positions x, y and z of the atoms in Angstrom, respectively.
In our example the first 67392 atom positions correspond to the atoms of subgroup 0 followed by the 33696 positions of subgroup 1 and the 16848 positions of subgroup 2.
**Wrong atom numbers specified for the subgroups and false ordering of the positions will result in false diffraction patterns.**

Exceeding the magnitude of 999.9999 Angstrom for x, y or z can decrease the accuracy, as only single precision numbers are used to store the atomic positions. 
For this reason, it is recommended to shift the cluster center to 0.

### Check your Parameters
In the _Windows Command Prompt_ (cmd) you can check if all parameters were read correct:

![Atomic Scattering Factors](/images/cmd_Para.PNG)

### Too Complicated, Try xyz2deb
More comprehensive and and easier to generate is a [xyz-file](/Data/UH3_7nm.xyz), where specific atomic ordering is no condition.  
Each line contains 6 arguments separated by whitspace, specifying the properties of a single atom/ion:

![xyz-file](/images/UH3xyz.PNG)

**The first argument** is the atom/ion type. This abbreviation should be consistent with the symbols used for the [atomic scattering factors](/Matlab/Tables/atomic_scattering_factors_RIETAN.txt).

**The second, third and fourth values** are the atom positons x, y and z, respectively. These values are given in Angstrom.

**The fourth and fifth values** are the occupancy and the isotropic temperature factor, respectively.

It is also possible to omit the isotropic temperature factors or both the isotropic temperature factors and the occupancies (see [*UH3_7nm_reduced.xyz*](Data/UH3_7nm_reduced.xyz)).
In such a case the occupancy is set to 1 and the isotropic temperature factor to 0.
This allows to use exported xyz-files from [_Vesta_](http://jp-minerals.org/vesta/en/) directly after the removal of the header. 
With an enhanced text editor providing column wise operations and enhanced replacement options, like [_Notepad++_](https://notepad-plus-plus.org/) occupancies and temperature facotors can be supplemented.

For the conversion of a xyz-file you can run the [_Matlab_](https://mathworks.com/) function [_xyz2deb.m_](/Matlab/xyz2deb.m):

![xyz2deb](/images/xyz2deb.PNG)

Alternatively, you can run the function via the command window without GUI by specifying the file, like *xyz2deb('C:\Users\rudolp2\Desktop\Cuda\cuDebye-v1.5\Data\UH3_7nm.xyz')*.
In the function [_xyz2deb.m_](/Matlab/xyz2deb.m): you can also edit the first lines, so your favourite properties are loaded automatically.
_You have not to worry about the maximum cluster size._

## Running cuDebye
First run the _Windows Command Prompt_ (cmd) and navigate to your folder containing _cuDebye.exe_, _cudart32_80.dll_ and _cudart64_80.dll_. Then you can start _cuDebye_ and pass a deb-file for the calculation:

- **cd** _PathWhereYouExtractedTheBranch_
- **cuDebye.exe** _PathOfDebFile_

![run](/images/run.PNG)

After the program has finished the calculation, you can check the output for warnings and errors (use fullt text search in some editor).

**It is recommended to keep the original folder structure from the branch**, but you can create new folders in the existing ones.
Nevertheless, you should pass on names of parent or children folders which containing folder names from the branch.

## Plotting
You can process and plot the simple output files with any program you like even with _Mircrosoft Excel_.
However, some [_Matlab_](https://mathworks.com/) functions are provided in the folder [_MatlabScripts_](/MatlabScripts).

### Plotting the Diffraction Patterns
