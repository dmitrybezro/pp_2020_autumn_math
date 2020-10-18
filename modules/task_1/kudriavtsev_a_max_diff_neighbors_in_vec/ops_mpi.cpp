// Copyright 2020 Kudriavtsev Alexander
#include <mpi.h>
#include <vector>
//#include <random>
#include <ctime>
#include <algorithm>
#include <iostream>
#include "ops_mpi.h"


std::vector<int> getRandomVector(int n) {
    std::vector<int> vec(n);
	srand((unsigned int)(time(0)));
    for (size_t  i = 0; i < n; ++i) vec[i] = rand() % 1000; 
    return vec;
}

int getSequentialOperations(std::vector<int>& vec) { // Is using & LEGAL? Probably yes, but Idk
    const size_t  n = vec.size();
    double t_b, t_e;//Time of the Beginning and of the End... (No End,No Beginning is a good song)
	t_b = MPI_Wtime();
	MyPair res;
	res.diff = abs(vec[1]-vec[0]);
	res.indx = 0;
	for(int i=1;(i+1)<n; ++i){
		int tmp = abs(vec[i+1]-vec[i]);
		if(tmp > res.diff) {
			res.diff = tmp;
			res.indx = i;
		}
	}
	t_e = MPI_Wtime();
	std::cout<<"Time of the Sequence :"<<t_e-t_b<<std::endl;
	
    return res.indx;
}

int getParallelOperations(std::vector<int>& global_vec) {
	//Beginning of COPYPASTA
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double t_b, t_e;
    if(rank == 0) t_b = MPI_Wtime();

    const int count_size_vector = global_vec.size();
    const int delta = count_size_vector / size;
    const int remed = count_size_vector % size;
    
    if (rank == 0) {
        for (int proc = 1; proc < size; ++proc) {
            if (proc <= remed) {
                MPI_Send(&global_vec[0] + proc * delta + proc - 2, // Is using &global_vec[0] LEGAL? CPPReference says YES since c++03 https://en.cppreference.com/w/cpp/container/vector
                    delta + 2, MPI_INT, proc, 0, MPI_COMM_WORLD); // MPI_Bcast (?) is faster
            }
            else {
                MPI_Send(&global_vec[0] + proc * delta + remed-1, delta+1, MPI_INT, proc, 0, MPI_COMM_WORLD);
            }
        }
    }

    std::vector<int> local_vec(delta);
    if (rank == 0) {
        local_vec = std::vector<int>(global_vec.begin(),
                                     global_vec.begin() + delta);
    } else {
        MPI_Status status;
        local_vec.resize(delta + size_t(rank > remed));
        MPI_Recv(&local_vec[0], delta + size_t(rank > remed), MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
	//End of COPYPASTA
    MyPair global_res;
    auto local_res = getSequentialOperations(local_vec); //Is this LEGAL? IDK, but 
    
    MPI_Reduce(&local_res, &global_res, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
    if (rank == 0) t_e = MPI_Wtime();
    std::cout << "Parallel time: " << t_e - t_b << std::endl;
    return global_res.diff;
}
