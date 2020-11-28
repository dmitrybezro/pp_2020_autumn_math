// Copyright 2020 Kudriavtsev Alexander
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <iostream>
#include "./ops_mpi.h"

std::vector<int> getRandomVector(int n) {  // Copied from https://github.com/allnes/pp_2020_autumn_math/pull/1
    std::vector<int> vec(n);
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    for (int i = 0; i < n; ++i) {
        vec[i] = gen() % 1000;
    }
    return vec;
}

int getSequentialOperations(std::vector<int> vec, int n) {  // Is using  & LEGAL? Travis says NO, but I say YES.
    double t_b, t_e;  // Time of the Beginning and of the End... (No End,No Beginning is a good song)
    t_b = MPI_Wtime();
    MyPair res;
    res.diff = abs(vec[1] - vec[0]);
    res.indx = 0;
    for (int i = 1; (i + 1) < n; ++i) {
        int tmp = abs(vec[i + 1] - vec[i]);
        if (tmp > res.diff) {
            res.diff = tmp;
            res.indx = i;
        }
    }
    t_e = MPI_Wtime();
    std::cout << "Time of the Sequence :" << t_e-t_b << std::endl;
    return res.indx;
}

int getParallelOperations(std::vector<int> global_vec, int count_size_vector) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double t_b, t_e;
    if (rank == 0) {
        t_b = MPI_Wtime();
    }

    const int delta = count_size_vector / size;
    const int remed = count_size_vector % size;
    if (rank == 0) {
        for (int proc = 1; proc < size; ++proc) {
            if (proc <= remed) {
                MPI_Send(&global_vec[0] + proc * delta + proc - 2,  // Is using &global_vec[0] LEGAL? YES
                    delta + 2, MPI_INT, proc, 0, MPI_COMM_WORLD);  // MPI_Bcast (?) is faster
            } else {  //  If an else has a brace on one side, it should have it on both... I don't like this.;
                MPI_Send(&global_vec[0] + proc * delta + remed - 1, delta + 1, MPI_INT, proc, 0, MPI_COMM_WORLD);
            }
        }
    }

    int* local_vec;
    const int n = (rank == 0) ? delta : (delta + 1 + static_cast<int>(rank <= remed));
    local_vec = new int[n];
    if (rank == 0) {
        for (int i = 0; i < delta; ++i) {
            local_vec[i] = global_vec[i];
        }
    } else {
        MPI_Status status;
        MPI_Recv(local_vec, n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
    MyPair global_res;
    MyPair local_res;

    local_res.diff = abs(local_vec[1] - local_vec[0]);
    local_res.indx = 0;
    for (int i = 1; (i + 1) < n; ++i) {
        int tmp = abs(local_vec[i + 1] - local_vec[i]);
        if (tmp > local_res.diff) {
            local_res.diff = tmp;
            local_res.indx = i;
        }
    }
    for (int i = 0; i < rank; ++i) {  // Not the best way
        local_res.indx += delta - static_cast<int>(i == 0) * 2 + static_cast<int>(i <= remed);
    }
    MPI_Reduce(&local_res, &global_res, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
    delete[] local_vec;
    if (rank == 0) {
        t_e = MPI_Wtime();
        std::cout << "Parallel time: " << t_e - t_b << std::endl;
    }
    return global_res.indx;
}
