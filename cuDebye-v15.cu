// cuDEBYE SOURCE CODE VERSION 1.5
// TO DO:
// - REWRITE TO DOUBLE PRECISION DISTANCE CALCULATIONS FOR BENCHMARKING
// - CONSIDER NOT CALLING SQRT (HISTOGRAM OF VALUE UNDER SQUARE -> problem with memory, no solution jet) IN KERNEL TO SAVE COMPUTATION TIME
// - USE INTEGER VALUES INSTEAD OF FLOAT AND CALCULATE IN FEMTO METERS INSTEAD OF ANGSTROM -> INTEGER OPERATIONS SHOULD REPLACE ROUND AND SINGLE PRECISION OPERATIONS WITH ACCEPTABLE ERROR
// - IMPLEMENT A CLEVER ALGORYTHM TO SET GRID AND BLOCK SIZE AUTOMATICALLY
// - BINARY FILE SUPPORT FOR FASTER INFORMATION EXCHANGE AND LESS MEMORY CONSUMPTION OR/AND PYTHON7MATLAB INTERFACE TO GET ARRAYS DIRECTLY
// - CREATE INTERFACE TO DISCUS (READ DISCUS STRUCTURES)
// - IMPLEMENT USAGE OF MORE GPU'S
// - MULTIPLE EMPTY LINES IN ASCII CAN CAUSE A CRASH DURING READING
// - HOST AND THRUST OPERATIONS ARE VERY INEFFICIENT (BUT FAST ENOUGH) -> MAYBE REWRITE THEM
// - ELIMINATE COMPILER WARNINGS FOR A MORE STABLE PROGRAM


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PREAMBLE: LIBARIES AND USEFULL BASIC FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Include cuda libaries for parallel computing
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include <sm_20_atomic_functions.h>

// Thrust libaries from the cuda toolkit for optimized vector operations
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/extrema.h>
#include <thrust/copy.h>
#include <thrust/gather.h>
#include <thrust/iterator/counting_iterator.h>

// Libaries for input and output streams for display results and read and write files.
// Better than the old printf shit
#include <fstream>		// File Stream
#include <iostream>		// Input/Output Stream
#include <iomanip>    // For Output Precision
#include <sstream>		// String Stream -> much easier than char[size] because size is set automatically
using namespace std;	// Normally all stream functions have to called via prefix std:: -> So functions can called withaout prefix (Example: std::cout -> cout)

						// Libary for measuring calculation time
#include <ctime>

						// define the mathematical constant pi
# define PI 3.14159265358979323846

						// Function to check if input file parsed via commandline exists
bool fexists(const char *filename) {
	ifstream file(filename);
	file.close();
	return bool(file.good());
}

// CUDA ERROR CHECK: can be wrapped around every device function, to abort and report if an error occurs.
// I just copied this from an guy of the Stack Overflow Community. 
// http://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ESSENTIAL KERNEL FUNCTIONS FOR FAST DISTANCE CALCULATION AND HISTOGRAM GENERATION
// INCLUDING A CONTROLLING HOST FUNCTION FOR PROPER KERNEL SELECTION
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// There is no hardware support for floating-point divisions (or integer divisions, for that matter) on the GPU,
// so these operations are implemented as software  subroutines that require additional registers for temporary
// storage.
// This is why 1/r is calculated by the host and passed as rp_rstep to the device and kernels -> saves a lot of computational time
// if thic calculation is avoided in the kernel.
// Additionally it should be noted that: 
// - The less operations within a kernel the faster the kernel is.
// - Single precision operations are faster than double precision operations and simple +,-,* operations are faster then complex math like power, sqrt, etc.
// - If a value is needed twice or more for some further calculation store this value in a local variable and access it to avoid a new calcualtion
//   in normal programs the time do not increase essentially but here this operation is called massivley multiple times which will result in increased
//   computing time.
// - Use pointers for larger arrys to prevent the transgfer of large arrays to the kernel each time it is called. But be aware of the fact that
//   the acess to the global memory is slower than to the local and shared memory of a kernel function.
// ...
// Kernel which will calculate and count all distances withaout checking if the passed atom indices exist.
__global__ void dist2(int *start, float *rp_rstep, float *x, float *y, float *z, unsigned int *nr) {
	// We pretend to have a part of distance matrix with atoms of the first subset in columns and the second subset in
	// the rows. Since all atom positions are stored together in three arrays, the kernel has to calculate the indcices
	// for the atom in row and column of the imaginary distance matrix. Therefore the start value depending on subgroup and
	// subset for the actual grid (see main and splitCalcDist function -- necessary because distance matrix can be larger than
	// grid -> see explanation above splitCalcDist) is necessary as well as the cuda built in variables
	// which return the block and grid indices to calculate the indices of the vector miming the matrix.
	// See: http://users.wfu.edu/choss/CUDA/docs/Lecture%205.pdf, Slide 16, "Flatten Matrices into linear Memory")
	int Idx_col = threadIdx.x + blockIdx.x*blockDim.x + start[0];
	int Idx_row = threadIdx.y + blockIdx.y*blockDim.y + start[1];
	// Calculating distances between the two atoms acessing the global memory to get thei coordinates.
	//float r = sqrt(pow(x[Idx_col] - x[Idx_row], 2) + pow(y[Idx_col] - y[Idx_row], 2) + pow(z[Idx_col] - z[Idx_row], 2));
	// The following shit is much faster
	float dx = x[Idx_col] - x[Idx_row];
	float dy = y[Idx_col] - y[Idx_row];
	float dz = z[Idx_col] - z[Idx_row];
	float r = dx*dx + dy*dy + dz*dz;
	// Calculate from the interatomic distance the channel/slot of nr which corresponds to the rounded distances value
	// with user defined accuracy/stepsize -> use r again to prevent a new variable definition.
	r = round(sqrt(r)*rp_rstep[0]);
	// Use atomic add inspired by fast histogram: https://devblogs.nvidia.com/parallelforall/gpu-pro-tip-fast-histograms-using-shared-atomics-maxwell/
	// At storage position of array nr corresponding to the calculated distance one is added -> for counting the frequency distances.
	// In case of simoultaneously the acces to one storeage place one theread has to wait until the other has finished.
	// Therefore the bottleneck here is the bandwith and the bus interface for global memory access, but in our case the time determining step
	// is the calculation, so this problem doesen't matter for the 980 GTX TI only an interface load of 5-10 percent is measured and the transfer
	// rate of 346 Gb/s guarantees a very short waiting time in the rare case of simultaniouse access and reduces the probaility of this case.
	// So the speedloss due to this is expected to be insignificant.
	atomicAdd(&nr[int(r)], 1);
}
// Kernel which do tha same as dist2, but checks for phantom atoms before calculation. In case of non defined atom index kernel returns without further
// calculation and is faster available for a new task.
__global__ void care_dist2(int *start, int *end, float *rp_rstep, float *x, float *y, float *z, unsigned int *nr) {
	int Idx_col = threadIdx.x + blockIdx.x*blockDim.x + start[0];
	int Idx_row = threadIdx.y + blockIdx.y*blockDim.y + start[1];
	if (Idx_col > end[0] || Idx_row > end[1]) {
		return;
	}
	float dx = x[Idx_col] - x[Idx_row];
	float dy = y[Idx_col] - y[Idx_row];
	float dz = z[Idx_col] - z[Idx_row];
	float r = dx*dx + dy*dy + dz*dz;
	r = round(sqrt(r)*rp_rstep[0]);
	atomicAdd(&nr[int(r)], 1);
}
// Kernel which do tha same as care_dist2 but for two equal atom subsets, so before the calculation it is checked that only the upper triangular
// distance matrix is calculated (col of subset larger than row of subset).
__global__ void dist(int *start, int *end, float *rp_rstep, float *x, float *y, float *z, unsigned int *nr) {
	int Idx_col = threadIdx.x + blockIdx.x*blockDim.x + start[0];
	int Idx_row = threadIdx.y + blockIdx.y*blockDim.y + start[1];
	if (Idx_col > Idx_row || Idx_row > end[1]) {
		return;
	}
	float dx = x[Idx_col] - x[Idx_row];
	float dy = y[Idx_col] - y[Idx_row];
	float dz = z[Idx_col] - z[Idx_row];
	float r = dx*dx + dy*dy + dz*dz;
	r = round(sqrt(r)*rp_rstep[0]);
	atomicAdd(&nr[int(r)], 1);
}

// Due to watchdog and to display the calculation process it is usefull to split the Grid into loops
// based on given grid and block size, to prevent to long kernel execution times. At a reasonable
// size speed loss do to loops is negligible
// The LoopX and LoopY mimes a grid for the grid, where each part isn't executed parallel but sucessively.
void splitCalcDist(int BlockDimXY, int GridDimXY, int *start_group, int *end_group, float *dev_rp_rstep, float *dev_x, float *dev_y, float *dev_z, unsigned int *dev_nr) {
	// For simplicity a quadratic block size is always used and a quadratic grid size is used to estimate
	// the necessary dimensions of the double loop and only changed for the last column and row loop.
	dim3 block(BlockDimXY, BlockDimXY);
	int GridDimX = GridDimXY;
	int GridDimY = GridDimXY;
	// Calculate number of atoms for subgroup a (columns) and b (rows)  
	int n_col = end_group[0] - start_group[0] + 1;
	int n_row = end_group[1] - start_group[1] + 1;

	// Calculate the loop dimension to cover all atoms of the column subgroup with quadratic grids.
	int LoopDimX = n_col / (GridDimXY*BlockDimXY); // integer division value will be floored
												   // If the number of atoms is a non multiple of GridDimXY*BlockDimXY not all distances will be covered. 
												   // Check this by calculating the remainder after devision
	int LastLoopGridDimX = n_col % (GridDimXY*BlockDimXY);
	if (LastLoopGridDimX != 0) {
		// If not every distance can be covered by a multiple of the defined grid size in X (col) direction,
		// another grid will be added to the X loop dimension and the integer LastLoopGridDimX will be changed
		// from 0 to the number of grids including the last incomplete grid (for a later if condition, so
		// the program knows when the quadratic grid size has to be modified to an rectangular size).
		// Additionally the grid dimension for the last column grid is calculated.
		LoopDimX++;
		LastLoopGridDimX = LastLoopGridDimX / BlockDimXY + 1;
	}
	// Do the same for row subgroup
	int LoopDimY = n_row / (GridDimXY*BlockDimXY);
	int LastLoopGridDimY = n_row % (GridDimXY*BlockDimXY);
	if (LastLoopGridDimY != 0) {
		LoopDimY++;
		LastLoopGridDimY = LastLoopGridDimY / BlockDimXY + 1; // integer devision floors the value so add 1 to ceil and have a grid covering all atoms
	}
	// Display the number of necessary loops to calculate all grids as well as the dimension of the last grid.
	cout << "\nNumber of Row Wise Loops: " << LoopDimY << ", Y-Grid-Dimension for the Last Loop: " << LastLoopGridDimY << endl;
	cout << "Number of Column Wise Loops: " << LoopDimX << ", X-Grid-Dimension for the Last Loop: " << LastLoopGridDimX << endl;
	cout << endl;

	// Introduce a boolean variable for the double loop so always a efficient kernel without
	// if conditions is called if this isn't necessary
	bool careX, careY;
	// Introduce host and device variables to storing the indices of atoms for each grid
	// -> this is necessary because the atom subgroups has to be splitted in subsets fullfilling
	// the conditions of the grid and block sizes
	int *start_subset = new int[2];
	int *end_subset = new int[2];
	int *dev_start_subset, *dev_end_subset;
	gpuErrchk(cudaMalloc((void **)&dev_start_subset, 2 * sizeof(int)));
	gpuErrchk(cudaMalloc((void **)&dev_end_subset, 2 * sizeof(int)));
	// Start loop over grids rowwise
	for (int LoopY = 0; LoopY < LoopDimY; LoopY++) {
		// Check if the current loop is the last rowwise loop for grid/subset calculation and set careY
		// to 1 and change the row grid dimension if it is the last loop.
		if (LoopY == LoopDimY - 1 && LastLoopGridDimY != 0) {
			GridDimY = LastLoopGridDimY;
			careY = 1;
		}
		else {
			GridDimY = GridDimXY;
			careY = 0;
		}
		// Start loop over columns for actual row
		for (int LoopX = 0; LoopX < LoopDimX; LoopX++) {
			// Check if the current loop is the last columnwise loop for grid/subset calculation and set careY
			// to 1 and change the column grid dimension if it is the last loop.
			if (LoopX == LoopDimX - 1 && LastLoopGridDimX != 0) {
				GridDimX = LastLoopGridDimX;
				careX = 1;
			}
			else {
				GridDimX = GridDimXY;
				careX = 0;
			}
			// Define the grid dimension for the actual row and column the double loop.
			dim3 grid(GridDimX, GridDimY);
			// Calculate the start and end indicies for the current column atom subset wich are parts of the current subgroups.
			start_subset[0] = start_group[0] + LoopX*GridDimXY*BlockDimXY;
			end_subset[0] = start_subset[0] + GridDimX*BlockDimXY - 1;
			// Check if calculated end of subset isn't larger than the end of subgroup -> necessary for last loop 
			// because grid can be a little to large due to ceil (see calculation of LastLoopGridDim).
			if (end_subset[0]>end_group[0]) { end_subset[0] = end_group[0]; }
			// Do the same  for the current row subset
			start_subset[1] = start_group[1] + LoopY*GridDimXY*BlockDimXY;
			end_subset[1] = start_subset[1] + GridDimY*BlockDimXY - 1;
			if (end_subset[1]>end_group[1]) { end_subset[1] = end_group[1]; }
			// Display current loop/grid (including size) and the Atom indices for which the distances are calculated
			cout << "Loop: " << LoopY << "," << LoopX << " -> Grid: " << grid.y << "x" << grid.x << " -> Atom Index: ";
			cout << start_subset[1] << "-" << end_subset[1] << ", ";
			cout << start_subset[0] << "-" << end_subset[0] << endl;
			// copy the atom indices of the current subset from the host to the device
			gpuErrchk(cudaMemcpy(dev_start_subset, start_subset, 2 * sizeof(int), cudaMemcpyHostToDevice));
			gpuErrchk(cudaMemcpy(dev_end_subset, end_subset, 2 * sizeof(int), cudaMemcpyHostToDevice));
			// Check some conditions to choose an efficient kernel function for distance calculation
			if (start_group[0] == start_group[1]) {
				// If the distances are calculated between the same subset (column atoms = row atoms) a kernel with permanent
				// if condition before calculation is launched to check that only the upper triangular distance matrix is
				// calculated. The if conditions which are executed for every kernel cost a lot of time but reduce the amount
				// of calculation to a half, which saves more time than necessary for the conditioning.
				dist <<<grid, block>>>(dev_start_subset, dev_end_subset, dev_rp_rstep, dev_x, dev_y, dev_z, dev_nr);
			}
			else if (careX == 0 && careY == 0) {
				// Standard calculation kernel for different subsets, where the row or column of the grid is quadratic as
				// defined by the user and no phantom atoms exist. This is usually the case if the last row or column loop/grid
				// isn't reached or even in case of the last loop if the subgroup was a multiple of block and grid size.
				// Becaus this is checked once before kernel launch, this kernel can run without multiple called if condition 
				// before distance calculation, which saves a lot of computational time since this operation has to be called for 80
				// 90 percent of disctances which are calculated. In case of several billion distances this can save a lot of time.
				// For smaller clusters the effect will be negligible.
				dist2 <<<grid, block>>>(dev_start_subset, dev_rp_rstep, dev_x, dev_y, dev_z, dev_nr);
			}
			else {
				// In case that the grid is larger than the actual atoms in row or column (is only used if a last loop
				// with non user defined grid size is necessary) an if conditioned kernel is laucnched which will test
				// before calculation if the atom truly exists or if it is an empty place (phantom atom) caused by the discrete increments
				// of the block size.
				care_dist2 <<<grid, block>>>(dev_start_subset, dev_end_subset, dev_rp_rstep, dev_x, dev_y, dev_z, dev_nr);
			}
			// Synchronize device to force the host to wait the kernel is finis´hed before a new is launched, otherwise we would get
			// an communication overhead.
			gpuErrchk(cudaDeviceSynchronize());
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTION TO REDUCE HISTOGRAM OF SUBGROUP SET AND CHECK FOR ERRORS IN CALCULATION
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Removes every entry of r and nr where nr=0 (squeeze histogram) and converts squeezed r and nr to double for 
// further calculations which require double precision (Debye double sum). Therfore thrust operations are used
// these are not as efficient as raw cuda but are easier to use and for the moment fast enough.
// ....
// Further informations according thrust see comments in thrust.cu
void reduceCheck(double *dev_r, unsigned int *dev_nr, int channels, int *start_group, int *end_group, unsigned int *dev_nrTmp, double *dev_r_reduced, double *dev_nr_reduced, int *channels_reduced, int row, int col, char *filename, int hist) {
	// Wrap begin of device arrays with thrust vector pointer trD_ indidcates the first storage place of the first value of the thrust vector,
	// therefore the type of stored values is necessary, so thrust knows how many bytes have to be read for one vector element.
	// tdR_+1 will point to the beginning byte of the second value.
	// Because these are only wraps modifing trD via thrust operations actually modifies the device arrays -> so basically trD_nr is the same as dev_nr.
	thrust::device_ptr<double> trD_r(dev_r);
	thrust::device_ptr<unsigned int> trD_nr(dev_nr);
	thrust::device_ptr<double> trD_r_reduced(dev_r_reduced);
	thrust::device_ptr<unsigned int> trD_nrTmp(dev_nrTmp);
	// Get maximum and minimum frequencies as well es their positions stored in a new thrust device vector.
	// Since the wrapped arrays contains no size information start (trD_nr) and endpoint (trD+number_of_elements) are necessary for vector operations.
	thrust::device_ptr<unsigned int> minimum = thrust::min_element(trD_nr, trD_nr + channels); // vector contains all minimum values (multiple times for each occurance)
	thrust::device_ptr<unsigned int> maximum = thrust::max_element(trD_nr, trD_nr + channels);
	// Get the value of minimum and maximum since they are all the same, just use first value.
	unsigned int min_value = minimum[0];
	unsigned int max_value = maximum[0];
	// Calculation of the position of the first occurance of minimum and maximum.
	int min_pos = minimum - trD_nr;
	int max_pos = maximum - trD_nr;
	// First check for errors with output of maximum and minimum values for user interpretation
	if (min_value < 0 || max_value < 0 || min_pos < 0 || min_pos >= channels || max_pos < 0 || max_pos >= channels) {
		cout << "WARNING: Something seems wrong with the following min or max values" << endl;
	}
	cout << "\nMinimum Frequency of Distances: " << min_value << " at " << min_pos << endl;
	cout << "Maximum Frequency of Distances: " << max_value << " at " << max_pos << endl;

	// Find indices of nonzero elements and copy the values to an reduced array
	// ...
	// Create vector to store index values for nonzero elements, which can be dynamically resized in contrary to a pointer vector.
	// But here we will use just a vector of the size of channels and do no resizing.
	// In contrary to a device pointer vector .begin() and .last() can be used to get storage place of first and last value. 
	// So basically the vector know its own length, it is smarter than a pointer vector, but passing it to function you have to pass the complete vector.
	thrust::device_vector<int> trD_Idx(channels);
	// Counting iterator from 1 (first) to channels (last) -> no vector but can be used as vector (creates values on the fly -> saves memory)
	thrust::counting_iterator<int> first(0);
	thrust::counting_iterator<int> last = first + channels;
	// Create an iterator called index iterator, this iterator can be used to count the copied values of an operation to another vector.
	typedef thrust::device_vector<int>::iterator IndexIterator;
	// Copy values from pseudo vector (first to last) in case that corresponding values of vector trD_nr are nonzero to trD_Idx
	// and store the storage place of the last copied number to vector trd_Idx as trD_Idx_end.
	IndexIterator trD_Idx_end = thrust::copy_if(first, last, // vector which is copied indicated by start and end
		trD_nr, // vector on which if condition is applied (should be the same length as the vector from which the values are copied)
		trD_Idx.begin(), // vector where the values are copied
		thrust::identity<unsigned int>()); // condition
										   // Calculate number of copied/nonzero elements which corresponds to unique distances and display them.
	channels_reduced[0] = trD_Idx_end - trD_Idx.begin();
	cout << "Number of Unique Distances: " << channels_reduced[0] << endl;
	// Copy values of r and nr for all indices of nonzero elements stored in trD_Idx to reduced device arrays. 
	// Actually they have the same size as r and nr reduced and every values beyond channels_reduced[0] is treated as data waste.  
	thrust::gather(trD_Idx.begin(), trD_Idx_end, trD_r, trD_r_reduced);
	thrust::gather(trD_Idx.begin(), trD_Idx_end, trD_nr, trD_nrTmp);
	// Copy the values of the reduced arrays containing the distance and the frequencies to the host.
	// Sum all frequencies to get the number of all distances from the subgroup set and convert nr_reduced from int to double
	// for further calculations (Debye double sum).
	double *r = new double[channels_reduced[0]];
	unsigned int *nr = new unsigned int[channels_reduced[0]];
	double *nr_reduced = new double[channels_reduced[0]];
	gpuErrchk(cudaMemcpy(r, dev_r_reduced, channels_reduced[0] * sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(nr, dev_nrTmp, channels_reduced[0] * sizeof(int), cudaMemcpyDeviceToHost));
	long long unsigned int nr_sum = 0; // for larger cluster number of all distances cause a integer overflow, therefore use int64
	for (int i = 0; i < channels_reduced[0]; i++) {
		nr_sum = nr_sum + long long unsigned int(nr[i]);
		nr_reduced[i] = double(nr[i]);
	}

	// Calculate and display number of unique values -> only for development process
	//unsigned int *nr_all = new unsigned int[channels];
	//gpuErrchk(cudaMemcpy(nr_all, dev_nr, channels * sizeof(int), cudaMemcpyDeviceToHost));
	//long long unsigned int nr_sum_all = 0;
	//int uniqueValues = 0;
	//for (int i = 0; i < channels; i++){
	//	if (nr_all[i] != 0){
	//		nr_sum_all = nr_sum_all + long long unsigned int(nr_all[i]);
	//		uniqueValues++;
	//	}
	//}

	// Simple error check for distance calculation.
	// ...
	// Calculate number of distances from subgroup set, and compare with nr_sum if thy don't match an integer overflow occured -> Cluster too large
	// for calculation or some other error. Same for the number of atoms and the observed zero distances.
	gpuErrchk(cudaMemcpy(dev_nr_reduced, nr_reduced, channels_reduced[0] * sizeof(double), cudaMemcpyHostToDevice));
	long long unsigned int n_atoms_col = long long unsigned int(end_group[0] - start_group[0] + 1);
	long long unsigned int n_atoms_row = long long unsigned int(end_group[1] - start_group[1] + 1);
	cout << "Number of Atoms in Row: " << n_atoms_row << endl;
	cout << "Number of Atoms in Column: " << n_atoms_col << endl;
	long long unsigned int nDistTh;
	int j, zeroDist;
	if (start_group[0] == start_group[1]) {
		nDistTh = (n_atoms_col*n_atoms_row - n_atoms_row) / 2 + n_atoms_col;
		if (nDistTh != nr_sum) {
			cout << "WARNING: " << nDistTh << " Distances Expected and " << nr_sum << " Observed!!!" << endl;
		}
		j = 0;
		zeroDist = 0;
		while (r[j] <= 0.0005) {
			zeroDist = zeroDist + nr[j];
			j++;
		}
		if (r[0] != 0 || nr[0] != n_atoms_col) {
			cout << "WARNING: " << "Number of Expected Zero Distances Missing or not Correct!!!" << endl;
			cout << "Number of Atoms within a Tolerance of 0.0005 A is " << zeroDist << endl;
		}
		while (r[j] < 0.25) {
			zeroDist = zeroDist + nr[j];
			j++;
		}
		if (zeroDist > n_atoms_col) {
			cout << "WARNING: " << "Possible Atom Collision!!! Found around " << zeroDist - n_atoms_col << " interatomic Distances lower than the Radius of an Hydrogen Atom." << endl;
		}
		cout << "Atoms Expected: " << n_atoms_col << " - Atoms Observed: " << nr[0] << endl;
	}
	else {
		nDistTh = n_atoms_col * n_atoms_row;
		if (nDistTh != nr_sum) {
			cout << "WARNING: " << nDistTh << " Distances Expected and " << nr_sum << " Observed!!!" << endl;
		}
		j = 0;
		zeroDist = 0;
		while (r[j] < 0.25) {
			zeroDist = zeroDist + nr[j];
			j++;
		}
		if (zeroDist > 0) {
			cout << "WARNING: " << "Possible Atom Collision!!! Found around " << zeroDist << " interatomic Distances lower than the Radius of an Hydrogen Atom." << endl;
		}
	}
	cout << "Distances Expected: " << nDistTh << " - Distances Observed: " << nr_sum << endl;


	// Writing the calculated diffraction data to an output file
	if (hist > 0) {
		clock_t time_start = clock();
		stringstream savefile; // convert the name of the input file to a stringstream
		savefile << filename << "x" << row << "x" << col; // append an I for 'Intensity' to this name to create the name of the output file
		cout << "\nWriting Distance File " << savefile.str() << endl;
		ofstream myfile; // Create an outputstrem named myfile
		myfile.open(savefile.str()); // use open property of outputstream to create/open the *.debI file -> this name is achieved by using the string property of the stringstream
		for (int i = 0; i < channels_reduced[0]; i++) { // write output to myfile instead to the cout (CommandWindow Out)
			myfile << setprecision(7) << r[i] << "\t" << nr[i] << endl;
		}
		myfile.close(); // close the file
						// Get end time for file reading
		clock_t time_end = clock();
		// Calculate and display time for reading the input file
		float elapsed_time = float(time_end - time_start) / CLOCKS_PER_SEC;
		cout << "Time for Writing Distance File: " << elapsed_time << " s" << endl;
	}


	// Synchronize device to force the host to wait the kernel is finis´hed before a new is launched, otherwise we would get
	// an communication overhead.
	gpuErrchk(cudaDeviceSynchronize());
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// KERNEL AND CONTROLLING FUNCTION TO CALCULATE SCATTERING PREFACTOR OF A SUBGROUP SET
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// I = ffoobb * sum(nr*sin(Kr)/Kr) .... simple sum over distances is the debye double sum over all atoms.
// Kerrnel function which gets atomic scattering factor, occupancy and isotropic thermal displacement (debye-waller) to calculate the factor
// ffoobb = fi*fj*Oi*Oj*exp[-Mi]*exp[-Mj] for every K/TwoTheta.
// M = 8*pi²*u²*sin²(theta)/lambda² = B*sin²(theta)/lambda², B = 8*pi²
__global__ void atomicScatter(int type1, int type2, int size_K, double *occ, double *beq, double *K, double *a, double *b, double *c, double *ffoobb) {
	// Kernel is executed for each K/TwoTheta (one dimensional grid)
	int Idx = blockIdx.x*blockDim.x + threadIdx.x;
	// Only execute if K/TwoTheta exists and is no phantom value, caused be discrete grid and block size.
	if (Idx < size_K) {
		double rp16pi2 = -0.006332573977646; // = (-1) * 1/(16*pi²)
		double negativeHalfSquaredS = K[Idx] * K[Idx] * rp16pi2; // = -sin²(theta)/lambda², s = 2*sin(theta)/lambda = 1/d
																 // Calculate occupancy and debye-waller part of the prefactor
		ffoobb[Idx] = occ[type1] * occ[type2];
		ffoobb[Idx] = ffoobb[Idx] * exp(negativeHalfSquaredS*(beq[type1] + beq[type2]));
		// Calculate atomic scattering factords from 11 parameter approximation.
		double f1 = c[type1];
		double f2 = c[type2];
		for (int i = 0; i < 5; i++) {
			f1 += a[type1 * 5 + i] * exp(b[type1 * 5 + i] * negativeHalfSquaredS);
			f2 += a[type2 * 5 + i] * exp(b[type2 * 5 + i] * negativeHalfSquaredS);
		}
		// Complement prefactor with calculated scattering factors
		ffoobb[Idx] = ffoobb[Idx] * f1*f2;
	}
}
// Controlling function to determine 1D block and grid size based one user defined 2D grid and block sizes
void atomicProp(int BlockDimXY, int type1, int type2, int nTT, double *dev_occ, double *dev_beq, double *dev_K, double *dev_a, double *dev_b, double *dev_c, double *dev_ffoobb) {
	dim3 block(BlockDimXY*BlockDimXY);
	dim3 grid(nTT / block.x + 1);
	cout << "\nCalculation of Atomic Properties" << endl;
	cout << "Grid: " << grid.y << "x" << grid.x << "\t";
	cout << "Block: " << block.y << "x" << block.x << endl;
	atomicScatter <<<grid, block>>>(type1, type2, nTT, dev_occ, dev_beq, dev_K, dev_a, dev_b, dev_c, dev_ffoobb);
	// Synchronize device to force the host to wait the kernel is finis´hed before a new is launched, otherwise we would get
	// an communication overhead.
	gpuErrchk(cudaDeviceSynchronize());
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// KERNEL AND CONROLLING FUNCTION FOR A ROWSUM OF AN IMAGINARY MATRIX
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Kernel to calculate the row sum of the following matrix to get the intensitys for each K/TwoTheta
//         r[0]    r[1]    r[2]  ....      =     Row Sum
// K[0]    I[0,0]  I[0,1]  I[0,2]				 Isum[0]
// K[1]    I[1,0]  I[1,1]  I[1,2]                Isum[1]
// K[3]    I[2,0]  I[2,1]  I[2,2]                Isum[2]
//   .
//   .
// ...
__global__ void debyeRowSum(int Kstart, int Kend, int size_nr, double *K, double *r, double *nr, double *ffoobb, double *Isum) {
	// Create an array to store double values in size of the block size (=BlockDimXY²) this value has to be set an given by the 
	// controlling function, if the size is define here as __shared__double partialSum[64] only a specific number and no variable can be used.
	extern __shared__ double partialSum[];
	// Calculate the one dimensional thread index as the block would be read column by column, because the shared array is only 1D
	int tid = threadIdx.y*blockDim.x + threadIdx.x;
	// Calculate row and column of the current thread, because the block is actually two dimensional. For the column Kstart isn't necessary, see
	// controlling function (grid is 1D as well as loop)
	int Idx_col = blockIdx.x*blockDim.x + threadIdx.x; // only for partial sum, these values will be added later
	int Idx_row = Kstart + blockIdx.y*blockDim.y + threadIdx.y; // row indicates which K/TwoTheta value is used
																// Check if K/TwoTheta exist or phantom value was generated by the discrete block and grid size.
	if (Idx_row <= Kend) {
		// Check if the first BlockDimX (col) values of interatomic distances are above 0.0005 A than it will be treated as distance between two different
		// atoms and if the value is below the program will assume it is the same atom and the deviation from zero was caused by a rounding error.
		// According to this calculate the value of the debye equation for the current K and r set and store it in the shared memory for the  current block.
		if (r[Idx_col] > 0.0005) {
			partialSum[tid] = 2 * nr[Idx_col] * sin(K[Idx_row] * r[Idx_col]) / (K[Idx_row] * r[Idx_col]); // times two because only distances of upper triangular distance matrix where calculated
		}
		else
		{
			partialSum[tid] = nr[Idx_col]; // sin(0)/0 = 1 so nr*sin(0)/0 = nr
		}
		// Since the grid dimension for columns is 1 only the values for r/nr up to BlockDimX would be calculated, therefore in this thread 
		// a for loop is started which adds on this inital value the block dimension until all distances are covered, which would have the
		// same threadIdx for further imaginary blocks in X (col) direction if the grid would be two dimensional. 
		// These values are just added to the partial sum of imaginary 1D thread tid. (see additional comments and schemes, RowSumSheme)
		Idx_col += blockDim.x;
		for (Idx_col; Idx_col < size_nr; Idx_col += blockDim.x) {
			partialSum[tid] += 2 * nr[Idx_col] * sin(K[Idx_row] * r[Idx_col]) / (K[Idx_row] * r[Idx_col]);
		}
	}
	// Wait until all threads have calculated the partial sum 
	__syncthreads();
	// For a sumation of rows only one thread is necessary for each row, so choose only threads with threadIdx.x=0 (first column in block).
	// Then simply sum the contents within a row of the block -> the result will be the complete row sum. 
	if (threadIdx.x == 0 && Idx_row <= Kend) {
		double Sum = partialSum[tid];
		for (int i = 1; i < blockDim.x; i++) {
			Sum += partialSum[tid + i];
		}
		// Calculate Ipart for current subgroup set
		Sum = Sum*ffoobb[Idx_row];
		// Add this value to the complet summed intensity for all atoms.
		Isum[Idx_row] += Sum;
	}
}

// Controlling function for row sum kernel to set block, grid and shared memory size and splits the calculastion/grid into 
// some loops to avoid to long kernel execution times (same procedure as in splitCalcDist but only one dimensional).
void intensity(int BlockDimXY, int GridDimXY, int size_K, int size_nr, double *dev_K, double *dev_r_reduced, double *dev_nr_reduced, double *dev_ffoobb, double *dev_Isum) {
	// Convert size integers in uint64 to prevent integer overflow during the calculation of LoopDimY.
	long long unsigned int uint64_size_K = long long unsigned int(size_K);
	long long unsigned int uint64_size_nr = long long unsigned int(size_nr);
	long long unsigned int uint64_BlockDimXY = long long unsigned int(BlockDimXY);
	long long unsigned int uint64_GridDimXY = long long unsigned int(GridDimXY);
	// Set the size of shared memory so that each thread of a block can store one double value in shared memory
	size_t sharedMem = BlockDimXY*BlockDimXY * sizeof(double);
	// Calculate the loop and grid dimension for rows, no dimension is set in X (col) direction, so that one block can do the complete row summation.
	// For this the user defined input size for an quadratic grid and block size is used.
	int LoopDimY = int((uint64_size_nr*uint64_size_K) / (uint64_BlockDimXY*uint64_BlockDimXY*uint64_GridDimXY*uint64_GridDimXY)) + 1;
	int GridDimY = (size_K / BlockDimXY + 1) / LoopDimY + 1;
	dim3 block(BlockDimXY, BlockDimXY); // block dimension is 2D to take advantage of shared memory for faster operation
	dim3 grid(1, GridDimY);
	// Show calculated loop and grid sizes
	cout << "\nEvaluation of the Debye Double Sum" << endl;
	cout << "Number of Unique Distances: " << size_nr << endl;
	cout << "LoopDimY: " << LoopDimY << " -> Grid: " << grid.y << "x" << grid.x << endl;
	// Introduce variables which store the start and end index for current loop/grid of K/TwoTheta 
	int start = 0;
	int end = 0;
	for (int LoopY = 0; LoopY < LoopDimY; LoopY++) {
		// Calculate start and end index of current grid and check if end index is larger than number of K/TwoTheta.
		start = LoopY*GridDimY*BlockDimXY;
		end = start + GridDimY*BlockDimXY - 1;
		if (end > size_K - 1) { end = size_K - 1; }
		// Display current loop and which rows/K are summed.
		cout << "Loop: " << LoopY << "\tK-Index: " << start << "-" << end << endl;
		// Launch kernel to calculate the matrix K[start:end] x r[all] and sum the rows. Beside grid and block size, the size of shared memory has
		// to be passed to the kernel as third value (for this a size_t type is necessary).
		debyeRowSum <<<grid, block, sharedMem>>>(start, end, size_nr, dev_K, dev_r_reduced, dev_nr_reduced, dev_ffoobb, dev_Isum);
		// Synchronize device to force the host to wait the kernel is finis´hed before a new is launched, otherwise we would get
		// an communication overhead.
		gpuErrchk(cudaDeviceSynchronize());
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BEGIN OF MAIN FUNCTION (READ FILE; STORE ATOM POSITIONS; ALLOCATE AND MEMORY ON GPU; COPY THE IMPORTANT DATA TO THE DEVICE; INITIATE CALCULATION)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Mostly pointers (indicated by *) are used to store values to avoid passing large arrays between functions
// and avoid problems with the device storage.
// Excample:
// int *E = new int[x] creates an 1D array with x storage places for integers (size: 4*x bytes)
// - The firstt four storage adresses contain the first value (the adresses of on pointer array are side by side)
// - Calling E will return the storage adress of the first value
// - Calling E[0] will return the first value
// - Calling &E[x-1] will return the first adress of the four storage adresses (4 bytes) of the last value in the pointer
//   -> the and will give the adresses even for normal declared variables
// - Overadressing E[x]=a or writing a on storage adress E+4*x will cause a manipulation of a other variable
//   or array depending on system reservation -> this could be a variable from a comlete other program or the 
//   System, which can cause a crash or of an own program variable, which will totally mess up any further 
//   calculation.

void main(int argc, char* argv[])
//void main()
{
	// Display Information
	cout << "\ncuDebye - Version: 1.5" << endl;
	cout << " - Calculation of distances with square root in the kernel function" << endl;
	cout << " - Writes intensity and distance files" << endl;
	cout << "-----------------------------------------------------------------------------" << endl;
	cout << "Author: Martin Rudolph" << endl;
	cout << "Technische Universit\x84t Bergakademie Freiberg" << endl;
	cout << "Institute of Materials Science" << endl;
	cout << "Gustav-Zeuner-Strasse 5" << endl;
	cout << "09599, Freiberg, Germany" << endl;
	cout << "E-Mail: m.s.rudolph@outlook.com\n" << endl;
	// Get FileName from the argument vector
	char *filename = argv[1];
	//char filename[] = "C:\\Users/rudolp2/Desktop/GammaAl2O3_10nm.deb";
	// Start timer for measuring time to load file
	clock_t time_start = clock();
	// Call function to check if input file exists (this function returns a boolean)
	if (fexists(filename) == 0) {
		cout << "\nERROR: >>" << filename << "<< does not exist" << endl;
		return;
	}
	// Define known variables stored in the first three lines of the input file
	double lambda, TwoThetaMin, TwoThetaMax, TwoThetaStep;
	double rmax, rstep;
	float rp_rstep;
	int CudaDevice, BlockDimXY, GridDimXY, hist, subgroups;
	// Create the reading stream "file" and initialize a string for a line by line scan
	ifstream file(filename);
	string line;

	// Before and in the header there should be no empty lines.

	// Get String of the first line and convert this string to an own temporary stream
	getline(file, line);
	istringstream tmp(line);
	// Extract the Index of the CUDA device which will be used for the calculation as well as the user
	// defined block and grid size in order to guard against the overstrain of an older GPU or to prevent 
	// watchdog crashes if the programm runs on windows and a display is connected to the selected 
	// calculating device or a kernel execution timeout is enabled.
	// If only one nVidia card is installed the device index is typically zero.
	// WARNING: THIS PROGRAMM CAN ONLY USE ONE GPU/DEVICE FOR CALCULATION; THE INDEX DOES NOT 
	// SPECIFY THE NUMBER OF PARALLEL COMPUTING GPUs; IT DEFINES THE EXECUTING GPU !!!
	tmp >> CudaDevice >> BlockDimXY >> GridDimXY;
	// Define the calculating cuda device specified by the input file
	gpuErrchk(cudaSetDevice(CudaDevice));
	// Initiate a structure to store the properties of an cuda device
	cudaDeviceProp prop;
	// Get the properties of the choosen calculating device
	gpuErrchk(cudaGetDeviceProperties(&prop, CudaDevice));
	// Because the free memory of the devices isn't stored in the properties structure initiate two additional
	// variables to get the free device memory.
	size_t mem_free = 0;
	size_t mem_tot = 0;
	gpuErrchk(cudaMemGetInfo(&mem_free, &mem_tot));
	// Display the important device properties
	cout << "\nCUDA PROPERTIES" << endl;
	cout << "----------------------------------------------------------------------------\n" << endl;
	cout << "Index of Selected Calculating Device: " << CudaDevice << endl;
	cout << "Number of Streaming Multiprocessors: " << prop.multiProcessorCount << endl;
	cout << "Clock Rate: " << prop.clockRate << " KHz" << endl;
	cout << "Kernel Timeout Enabled: " << prop.kernelExecTimeoutEnabled << endl;
	// Thread and block properties of the device
	cout << "\nMax Threads Per Block: " << prop.maxThreadsPerBlock << endl;
	cout << "Max Threads per Block in X-Dimension: " << prop.maxThreadsDim[0] << endl;
	cout << "Max Threads per Block in Y-Dimension: " << prop.maxThreadsDim[1] << endl;
	cout << "Threads per Warp: " << prop.warpSize << endl; // a block is actually splited into warps, therfore the size of an block should be ideally a multiple of 32
	cout << "\nMax Blocks per Grid in X-Dimension: " << prop.maxGridSize[0] << endl;
	cout << "Max Blocks per Grid in Y-Dimension: " << prop.maxGridSize[1] << endl;
	cout << "\nNote: in Cuda the X describes the Column and the Y the Row. ";
	cout << "As it is typical in C, Python and Matlab for displaying the Grid and Block Sizes ";
	cout << "the following format is used: Rows x Columns." << endl;
	// Available memory
	cout << "\nTotal Global Memory: " << prop.totalGlobalMem << " byte" << endl;
	cout << "Free Available Memory: " << mem_free << " byte" << endl;
	cout << "Shared Memory per Block: " << prop.sharedMemPerBlock << " byte" << endl;
	// User defined block and grid size
	cout << "\nUSER DEFINED GRID AND BLOCK SIZE" << endl;
	cout << "----------------------------------------------------------------------------\n" << endl;
	cout << "Block: " << BlockDimXY << "x" << BlockDimXY << "(" << BlockDimXY*BlockDimXY << ")";
	cout << "\tGrid: " << GridDimXY << "x" << GridDimXY << "(" << GridDimXY*GridDimXY << ")" << endl;

	// Get string from the second line, clear the temporary string stream and load the new string into this
	// stream. Extract from this stream the wavelength the lower and upper limit for 2Theta as well as the
	// step to specify for which range the diffractogram  is calculated.
	getline(file, line);
	tmp.clear();
	tmp.str(line);
	tmp >> lambda >> TwoThetaMin >> TwoThetaMax >> TwoThetaStep;
	int nTT = int((TwoThetaMax - TwoThetaMin) / TwoThetaStep) + 1; // number of TwoTheta discrete TwoTheta values
																   // Create pointers to store TwoTheta and K at the Host
	double *TwoTheta = new double[nTT];
	double *K = new double[nTT];
	// Calculate TwoTheta as well as K from the given input values
	for (int i = 0; i < nTT; i++) {
		TwoTheta[i] = TwoThetaMin + double(i)*TwoThetaStep;
		K[i] = 4 * PI*sin(TwoTheta[i] / 2 * PI / 180) / lambda;
	}
	// Display settings for calculation of the diffractogram
	cout << "\nDIFFRACTOGRAM SETTINGS" << endl;
	cout << "----------------------------------------------------------------------------\n" << endl;
	cout << "Wavelength: " << lambda << " A" << endl;
	cout << "2ThetaMin: " << TwoThetaMin << "\370" << "\t2ThetaMax: " << TwoThetaMax << "\370" << "\t2ThetaStep: " << TwoThetaStep << "\370" << endl;
	cout << "Number of 2Theta values: " << nTT << "\tStart: 0 = " << TwoTheta[0] << "\370" << "\tEnd: " << nTT - 1 << " = " << TwoTheta[nTT - 1] << "\370" << endl;
	cout << "k = 4*pi*sin(theta)/lambda: " << "\tStart: " << K[0] << "\tEnd: " << K[nTT - 1] << endl;

	// Read from the third line the maximum possible distance between two atoms within the cluster
	// as well as the accuracy for the calculation of interatomic distances (in Angstrom). 
	// Usually the rmax value should be a bit larger than the real maximum to prevent a crash in case of 
	// rounding errors. I usually add 5% and round this number up to the next integer. 
	// For the accuracy the perfect value seems to be 0.0001. Be Carefull above 0.001 the diffractogram
	// will get worse and below 0.0001 you will see no real effect due to the single precission calculation
	// and the fact, that atomic position in structure files are in most cases limited to 4 digits due to uncertainties.
	// Anyway it should be noted, that due to the single precision operations to speed up the calculation the
	// efficient accuracy is not 0.0001 it should be around 0.0003 as long as the maximum distance between 
	// two atoms is less than 999 A -> At the time I haven't checked this properly.
	// The issue is:
	// - Single precission allows 7 significant nearly exact digits -> 999.1234
	// - Error = (x*dx+y*dy+z*dz)/sqrt(x^2+y^2+z^2)
	// - Assume x=y=z and dx=dy=dz=0.00005 -> Error sqrt(3)*dx = 0.000087
	// - Rounding Error due to 7-8 significant digits neglecteted for the square sum x^2+y^2+z^2
	//   -> therefore a value around 0.0003 is assumed as worst case error
	// - Clusters above 1000 A can cause larger errors because dx_max rises to nearly 0.0005 due to the
	//   precision of single values as well as the rounding error of the square sum becomes more problematic
	// Please look carefully at your diffractogram for larger clusters.
	getline(file, line);
	tmp.clear();
	tmp.str(line);
	// Extract from the third line string/stream rmax and rstep
	tmp >> rmax >> rstep >> hist;
	// Calculate r and the number of cahnnels/bins to store the frequencies of interatomic distances
	int channels = int(rmax / rstep) + 2;
	double *r = new double[channels];
	for (int i = 0; i < channels; i++) {
		r[i] = double(i)*rstep;
	}
	// Define reciprocal rstep for global parallel device functions because a division within such a 
	// function is problematic and time consuming.
	rp_rstep = float(1 / rstep);
	// Display settings for interatomic distance histogram
	cout << "\nDISTANCE CALCULATION SETTINGS" << endl;
	cout << "----------------------------------------------------------------------------\n" << endl;
	cout << "Rmax: " << rmax << " A\tRstep: " << rstep << " A" << endl;
	cout << "Number of Channels: " << channels << "\tStart: 0 = " << r[0] << " A\tEnd: " << channels - 1 << " = " << r[channels - 1] << " A" << endl;
	if (hist == 1) {
		cout << "Output Distance File: True \t Output Intensity File: False" << endl;
	}
	else if (hist == 2) {
		cout << "Output Distance File: True \t Output Intensity File: True" << endl;
	}
	else {
		hist = 0;
		cout << "Output Distance File: False \t Output Intensity File: True" << endl;
	}

	// Extract the number from the fourth line in file which should contain the number of atomic subgroups of 
	// identical type, occupancy and isotropic temperature factor
	getline(file, line);
	tmp.clear();
	tmp.str(line);
	tmp >> subgroups;
	cout << "\nCLUSTER AND ATOMIC INFORMATIONS" << endl;
	cout << "----------------------------------------------------------------------------\n" << endl;
	cout << "Number of Subgroups: " << subgroups << endl;
	// Based on the number of subgroups the next lines are read. Each line should contain the number of the atoms,
	// the atom type, the occupancy and the isotropic temperature factor and the 11 parameter atomic scattering factor
	// approximation in the subgroup. Beside the extraction of
	// these informations. The indices of atoms for beginning and ending of subgroups are calculateted. This is 
	// necessary, because all the x,y,z-positions will be stored in three arrays for all subgroups together.
	int *n_atoms = new int[subgroups]; // array which stores the atom number of each subgroup
	int *start = new int[subgroups]; // array which stores the atom index where a subgroup starts
	int *end = new int[subgroups]; // array which stores the atom index where a subgroup starts
	string *type = new string[subgroups];
	double *occ = new double[subgroups];
	double *beq = new double[subgroups];
	double *a = new double[subgroups * 5]; // arrays to store scattering factors for all subgroups
	double *b = new double[subgroups * 5];
	double *c = new double[subgroups];
	start[0] = 0; // first subgroup starts at 0
				  // Extract information from every subgroup and calculate the start and end indices of each subgroup and
				  // display them for user control.
	for (int i = 0; i < subgroups; i++) {
		getline(file, line);
		tmp.clear();
		tmp.str(line);
		tmp >> n_atoms[i] >> type[i] >> occ[i] >> beq[i] >> a[i * 5 + 0] >> a[i * 5 + 1] >> a[i * 5 + 2] >> a[i * 5 + 3] >> a[i * 5 + 4] >> b[i * 5 + 0] >> b[i * 5 + 1] >> b[i * 5 + 2] >> b[i * 5 + 3] >> b[i * 5 + 4] >> c[i];
		cout << "\nGroup " << i << " - Number of Atoms: " << n_atoms[i] << endl;
		end[i] = start[i] + n_atoms[i] - 1;
		if (i < subgroups - 1) {
			start[i + 1] = start[i] + n_atoms[i];
		}
		cout << "Atom index from " << start[i] << " to " << end[i] << endl;
		cout << "Atom Type: " << type[i] << ", Occupancy: " << occ[i] << ", Beq: " << beq[i] << endl;
		cout << "Scattering Coefficients:" << endl;
		cout << "a[1-5]: ";
		for (int j = 0; j < 5; j++) {
			cout << a[i * 5 + j] << ", ";
		}
		cout << "\nb[1-5]: ";
		for (int j = 0; j < 5; j++) {
			cout << b[i * 5 + j] << ", ";
		}
		cout << "\nc: " << c[i] << endl;
	}
	int n_atoms_sum = end[subgroups - 1] + 1; // calculate number of all atoms

											  // Based on the total number of atoms three 1D arrays using pointers (*) are created to store the atomic 
											  // positions unsing single precision.
	float *x = new float[n_atoms_sum];
	float *y = new float[n_atoms_sum];
	float *z = new float[n_atoms_sum];
	cout << "\nMemory allocated for " << n_atoms_sum << " Atoms." << endl;
	cout << "Reading Atomic Data..." << endl;
	// Read line by line x, y and z position of the atoms and ignore empty lines, which can be added
	// in the input file to show clearly start and en of a subgroup (not necessary).
	// Be carefull empty lines followed by empty lines as well as to much empty lines at the end
	// of the file can cause a program crash -> Why is not clear at the moment.
	int i = 0;
	while (getline(file, line))
	{
		if (line.empty() == 0) // if line is non empty read atomic positions and stor them
		{
			istringstream tmp(line);
			tmp >> x[i] >> y[i] >> z[i];
			i++;
		}
	}
	// Show number of loaded atoms and check if they match the number specified by the subgroups in the header.
	cout << i << " Atoms are loaded!" << endl;
	if (i != n_atoms_sum) {
		cout << "\n\nERROR: Number of extracted Atoms doesn't match the number specified in the preamble!" << endl;
		return;
	}
	// close file
	file.close();
	// Get end time for file reading
	clock_t time_end = clock();
	// Calculate and display time for reading the input file
	float elapsed_time = float(time_end - time_start) / CLOCKS_PER_SEC;
	cout << "\nTime for Loading the Input File: " << elapsed_time << " s" << endl;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Allocate and Intialize important Variables on the GPU
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start clock for measuring allocation time
	time_start = clock();
	// Get and show free device memory before memory allocation on the GPU
	cout << "\nALLOCATING THE NECESSARY MEMORY ON THE DEVICE" << endl;
	cout << "----------------------------------------------------------------------------\n" << endl;
	gpuErrchk(cudaMemGetInfo(&mem_free, &mem_tot));;
	cout << "Free Memory before Allocation: " << mem_free << " byte" << endl;

	// Initialize pointer adresses stored in the host but pointing to the device memory
	float *dev_rp_rstep, *dev_x, *dev_y, *dev_z;
	double *dev_r_reduced, *dev_nr_reduced, *dev_r, *dev_K, *dev_occ, *dev_beq, *dev_a, *dev_b, *dev_c, *dev_ffoobb, *dev_Isum;
	unsigned int *dev_nr, *dev_nrTmp;
	int *channels_reduced = new int[1]; // pointer for single varible to store the number of reduced channels (see function reduceCheck)

										// Reserve storage place in bytes on the device, with no specific storage type (void**) as typical for the host
										// Therefor the pointer from the host is used and to provide enougth bytes to store different types 
										// additional the number of storage places/bytes have to be passed. This can easiely calculateted from
										// by length_of_array*size_of_type
										// ...
										// Reciprocal accuaracy of distances, float type -> see explaination at atom position allocation on device.  
	gpuErrchk(cudaMalloc((void**)&dev_rp_rstep, sizeof(float))); // only one float
																 // Array of discrete interatomic diatances and array to stor their frequencies.
	gpuErrchk(cudaMalloc((void**)&dev_r, channels * sizeof(double))); // number of discrete r * 8 byte
																	  // nr as integer for histogram/frequency, since only the number of occuring distances is counted
																	  // via atomic add unsigned integer is necessary for maximum performance. This allows a counting up
																	  // to 2^32 unique distances which should be enougth even for larger clusters.
																	  // For temporary operations an additional variable is created.
	gpuErrchk(cudaMalloc((void**)&dev_nr, channels * sizeof(int)));
	gpuErrchk(cudaMalloc((void**)&dev_nrTmp, channels * sizeof(int)));
	// Creating reduced variables to store only values where nr!=0.
	// Use same length (channels) as nr for worst case that every distance r occurs at least once.
	// They are stored from index 0 to channels_reduced[0]
	// nr contains only the distances from atom i to j and has to be multiplied with 2 to mimic the distance from j to i.
	// For Calculation of Intensity double precision is necessary (see intensity and debyeRowSum), also in case
	// that nr for one distance is larger than 1/2*2^32 this format is necessary to prevent integer overflow
	// when multiplying with 2.
	gpuErrchk(cudaMalloc((void**)&dev_r_reduced, channels * sizeof(double)));
	gpuErrchk(cudaMalloc((void**)&dev_nr_reduced, channels * sizeof(double)));
	// K = 4*pi*sin(Theta)/lambds and array to store corresponding intensity
	gpuErrchk(cudaMalloc((void**)&dev_K, nTT * sizeof(double))); // nTT -> number of TwoTheta bzw. K values
	gpuErrchk(cudaMalloc((void**)&dev_Isum, nTT * sizeof(double)));
	// Atom positions float precision for faster calculation of distance in kernel functions.
	gpuErrchk(cudaMalloc((void**)&dev_x, n_atoms_sum * sizeof(float)));
	gpuErrchk(cudaMalloc((void**)&dev_z, n_atoms_sum * sizeof(float)));
	gpuErrchk(cudaMalloc((void**)&dev_y, n_atoms_sum * sizeof(float)));
	// Occupancies and Beq and atomic scattering factors only necessary for intensity calculation, 
	// therefore double precision necessary.
	gpuErrchk(cudaMalloc((void**)&dev_occ, subgroups * sizeof(double)));
	gpuErrchk(cudaMalloc((void**)&dev_beq, subgroups * sizeof(double)));
	// Atomic scattering factors
	gpuErrchk(cudaMalloc((void**)&dev_a, 5 * subgroups * sizeof(double))); // *5 beacause every subgroup has 5 a parameters
	gpuErrchk(cudaMalloc((void**)&dev_b, 5 * subgroups * sizeof(double)));
	gpuErrchk(cudaMalloc((void**)&dev_c, subgroups * sizeof(double)));
	// Array to store fi*fj*Occi*Occj*DebWalli*DebWallj for every TwoTheta -> see atomicProp and atomicScatter 
	gpuErrchk(cudaMalloc((void**)&dev_ffoobb, nTT * sizeof(double)));

	// Copy host values to the device and set counting values to zero
	gpuErrchk(cudaMemcpy(dev_rp_rstep, &rp_rstep, sizeof(float), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(dev_r, r, channels * sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(dev_K, K, nTT * sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemset(dev_nr, 0, channels * sizeof(int)));
	gpuErrchk(cudaMemset(dev_Isum, 0, nTT * sizeof(double)));
	gpuErrchk(cudaMemcpy(dev_x, x, n_atoms_sum * sizeof(float), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(dev_y, y, n_atoms_sum * sizeof(float), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(dev_z, z, n_atoms_sum * sizeof(float), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(dev_occ, occ, subgroups * sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(dev_beq, beq, subgroups * sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(dev_a, a, 5 * subgroups * sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(dev_b, b, 5 * subgroups * sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(dev_c, c, subgroups * sizeof(double), cudaMemcpyHostToDevice));

	// Get free memory of device after allocation and display it.
	gpuErrchk(cudaMemGetInfo(&mem_free, &mem_tot));
	cout << "Free Memory after Allocation: " << mem_free << " byte" << endl;
	// Calculate and display time for allocation
	time_end = clock();
	elapsed_time = float(time_end - time_start) / CLOCKS_PER_SEC;
	cout << "\nTime for Allocating the necessary Memory on the Device: " << elapsed_time << " s" << endl;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Calculation of each Subgroup Combination (Distances+Intensity)
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start timer for actual calculation
	time_start = clock();
	cout << "\n\n----------------------------------------------------------------------------" << endl;
	cout << "START OF ACTUAL CALCULATION" << endl;
	cout << "----------------------------------------------------------------------------" << endl;
	// Define variables to store start and end of atom indices for current subgroup combination. 
	int *current_group_start = new int[2];
	int *current_group_end = new int[2];
	// Run double loop over each possible subroup combination, imitating a matrix via loops
	//     0    1    2    3
	// 0  00   01   02   03
	// 1       11   12   13
	// 2            22   23
	// 3                 33
	for (int row = 0; row < subgroups; row++) {
		for (int col = row; col < subgroups; col++) {
			// Get start and end atom index for subgroup combination of the pseudo matrix 
			current_group_start[0] = start[col]; // col in cuda X
			current_group_start[1] = start[row]; // row in cuda Y
			current_group_end[0] = end[col];
			current_group_end[1] = end[row];
			// Display the current subgroup for Calculation
			cout << "\n\n>>>> Calculation of Subgroup Number " << row << " and " << col << endl;
			// Calculate the histogram of distances -> most time consuming process
			splitCalcDist(BlockDimXY, GridDimXY, current_group_start, current_group_end, dev_rp_rstep, dev_x, dev_y, dev_z, dev_nr);
			// Reduce the arrays r and nr to r_reduced and nr_reduced in combination with an error check for 
			// the calculation. (could be improved but it is fast enough compared to the splitCalc and intensity) 
			reduceCheck(dev_r, dev_nr, channels, current_group_start, current_group_end, dev_nrTmp, dev_r_reduced, dev_nr_reduced, channels_reduced, row, col, filename, hist);
			if (hist != 1) {
				// Calculate the preamble of each subgroup for TwoTheta containign the atomic scattering factors,
				// the occupancies as well as the isotropic temperature factors. (could be improved but it is fast enough compared to the splitCalc and intensity)
				atomicProp(BlockDimXY, row, col, nTT, dev_occ, dev_beq, dev_K, dev_a, dev_b, dev_c, dev_ffoobb);
				// Calculate the portiton of scattered intensity from the actual subgroup combinations, and add
				// these part to the intensity of the previouse parts -> second most time consuming operation,
				// which is small compared to splitCalc but increse as the number of unique distances increases.
				intensity(BlockDimXY, GridDimXY, nTT, channels_reduced[0], dev_K, dev_r_reduced, dev_nr_reduced, dev_ffoobb, dev_Isum);
				// Reset histogram counting array for next subgroup combination.
			}
			gpuErrchk(cudaMemset(dev_nr, 0, channels * sizeof(int)));
		}
	}
	// Calculate and display the time for the time intense calculation stuff
	time_end = clock();
	elapsed_time = float(time_end - time_start) / CLOCKS_PER_SEC;
	cout << "\nTime for Calculation of Distances and\nEvaluation of the Debye Double Sum: " << elapsed_time << " s" << endl;
	if (hist != 1) {
		// Create host variable and transfer the calculated intensity from the device to the host.
		double *Isum = new double[nTT];
		gpuErrchk(cudaMemcpy(Isum, dev_Isum, nTT * sizeof(double), cudaMemcpyDeviceToHost));
		// Display the first 10 and last 10 values of the calculated diffraction data.
		cout << "\n\n" << endl;
		cout << "2Theta\tIntensity" << endl;
		cout << "------------------" << endl;
		for (int i = 0; i < 10; i++) {
			cout << TwoTheta[i] << "\t" << Isum[i] << endl;
		}
		cout << "\n...\n" << endl;
		for (int i = nTT - 11; i < nTT; i++) {
			cout << TwoTheta[i] << "\t" << Isum[i] << endl;
		}

		// Writing the calculated diffraction data to an output file
		time_start = clock();
		cout << "\n Writing Output File..." << endl;
		stringstream savefile; // convert the name of the input file to a stringstream
		savefile << filename << "I"; // append an I for 'Intensity' to this name to create the name of the output file
		ofstream myfile; // Create an outputstrem named myfile
		myfile.open(savefile.str()); // use open property of outputstream to create/open the *.debI file -> this name is achieved by using the string property of the stringstream
		for (int i = 0; i < nTT; i++) { // write output to myfile instead to the cout (CommandWindow Out)
			myfile << TwoTheta[i] << "\t" << Isum[i] << endl;
		}
		myfile.close(); // close the file
		time_end = clock();
		elapsed_time = float(time_end - time_start) / CLOCKS_PER_SEC;
		cout << "Time for Writing Output File: " << elapsed_time << " s" << endl;
	}
	return;
}