# Parallel Image Processing using MPI (C, Linux)

This project demonstrates parallel image processing using MPI (Message Passing Interface) in C.

## Context
This project was developed as part of a university laboratory for Parallel Systems Architecture.  
The implementation follows concepts and guidance provided during lab sessions.

## Features
- Parallel processing of PGM images
- Domain decomposition across multiple processes
- 2D Cartesian topology (MPI_Cart_create)
- Inter-process communication using:
  - MPI_Send
  - MPI_Recv
  - MPI_Sendrecv
- Halo exchange between neighboring processes
- Custom MPI datatypes (MPI_Type_vector)

## Technologies
- C programming language
- MPI (e.g. OpenMPI / MPICH)
- Linux (Ubuntu)

## How it works
The image is divided into blocks distributed across multiple processes.  
Each process computes its portion and exchanges boundary data (halo regions) with neighbors to ensure correct results.

## Compile
```bash
mpicc main.c -o app
```

## Run
```bash
mpirun -np 4 ./app
```

## Input/Output
- Input: PGM image file
- Output: processed PGM image

## Notes
- Implemented during laboratory sessions (not a fully independent project)
- Focused on understanding MPI concepts and parallel data distribution
