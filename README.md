# Multi-Objective Shortest Path Algorithms

This project implements static and dynamic multi-objective shortest path (MOSP) algorithms in both serial and parallel (MPI + OpenMP) C++. It also includes a Python script for generating large synthetic graph datasets.

## Project Structure

```
Source/
  Serial/
    main.cpp        # Serial C++ implementation
    graph.txt       # Example graph data
  Parallel/
    main.cpp        # Parallel C++ implementation (MPI + OpenMP)
    graph.txt       # Example graph data
Data/
  generate_dataset.py # Python script to generate datasets
  data1.txt           # Example large dataset
```

## Features
- **Serial and Parallel Implementations:**
  - Serial version in `Source/Serial/main.cpp`.
  - Parallel version in `Source/Parallel/main.cpp` using MPI and OpenMP.
- **Static and Dynamic MOSP:**
  - Handles both static and dynamic changes (edge insertions/removals) in the graph.
- **Flexible Input:**
  - Reads graph data from text files.
- **Dataset Generation:**
  - Python script to generate large random graph datasets.

## Input Format
Each line in a graph file represents an edge:
```
source_node,destination_node,weight1,weight2,...,weightN,
```
- The file ends with a `~` character.
- Example (2 objectives):
  ```
  0,1,4,5,
  1,2,3,7,
  ...
  4,6,7,8,~
  ```

## How to Build and Run

### Serial Version
1. **Navigate to the serial source directory:**
   ```sh
   cd Source/Serial
   ```
2. **Compile:**
   ```sh
   g++ -std=c++11 -o mosp_serial main.cpp
   ```
3. **Run:**
   ```sh
   ./mosp_serial
   ```
   - Uses `graph.txt` in the same directory by default.

### Parallel Version (MPI + OpenMP)
1. **Navigate to the parallel source directory:**
   ```sh
   cd Source/Parallel
   ```
2. **Compile:**
   ```sh
   mpicxx -fopenmp -std=c++11 -o mosp_parallel main.cpp
   ```
3. **Run (example with 4 processes):**
   ```sh
   mpirun -np 4 ./mosp_parallel
   ```
   - Uses `graph.txt` in the same directory by default.

### Dataset Generation
1. **Navigate to the Data directory:**
   ```sh
   cd Data
   ```
2. **Run the Python script:**
   ```sh
   python generate_dataset.py
   ```
   - This will generate a file named `data.txt` with 10,000,000 edges and 10,000 nodes by default.
   - You can modify the script to change the number of rows, nodes, or weight range.

## Requirements
- **Serial Version:** C++11 compatible compiler
- **Parallel Version:** MPI (e.g., OpenMPI), OpenMP, C++11 compatible compiler
- **Dataset Generation:** Python 3.x

## Notes
- The number of objectives (weights per edge) is set in the source code (`number_of_objectives`). Ensure your input files match this number.
- The example `graph.txt` files are small; use the Python script to generate larger datasets for benchmarking.
- The code demonstrates both static and dynamic MOSP, including edge insertions and removals.

## Author
Faris Ali
Faiz ul Hassan
Jawad Ahmad Khan