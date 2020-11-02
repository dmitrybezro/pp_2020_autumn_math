// Copyright 2020 Lebedev Andrew
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_1/lebedev_a_matrix_max_elem/matrix_max_elem.h"


std::vector<int> getRandomMatrix(int mtx_size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> mtx(mtx_size);
    for (int  i = 0; i < mtx_size; i++) {
            mtx[i] = gen() % 100;
    }
    return mtx;
}

int getSequentialOperations(std::vector<int> mtx) {
    const int  size = mtx.size();
    if (size == 0) {
        return -1;
    }
    int max_elem = mtx[0];
    for (int i = 1; i < size; i++) {
        max_elem = std::max(max_elem, mtx[i]);
    }
    return max_elem;
}

int getParallelOperations(std::vector<int> global_mtx, int mtx_size) {
    double time_start, time_stop;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        time_start = MPI_Wtime();
    }
    if (mtx_size == 0) {
        return -1;
    }
    const int delta = mtx_size / size;
    const int rem = mtx_size % size;
    int global_max;
    if (rank == 0) {
        for (int process = 1; process < size; process++) {
            if (process <= rem) {
                MPI_Send(&global_mtx[0] + process * delta + process - 1,
                    delta + 1, MPI_INT, process, 0, MPI_COMM_WORLD);
            } else {
                MPI_Send(&global_mtx[0] + process * delta + rem,
                    delta, MPI_INT, process, 0, MPI_COMM_WORLD);
            }
        }
    }
    std::vector<int> local_mtx(delta);
    if (rank == 0) {
        local_mtx = std::vector< int>(global_mtx.begin(),
            global_mtx.begin() + delta);
    } else {
        MPI_Status status;
        if (rank <= rem) {
            local_mtx.resize(delta + 1);
            MPI_Recv(&local_mtx[0], delta + 1, MPI_INT,
                0, 0, MPI_COMM_WORLD, &status);
        } else {
            MPI_Recv(&local_mtx[0], delta, MPI_INT,
                0, 0, MPI_COMM_WORLD, &status);
        }
    }
    int local_max = getSequentialOperations(local_mtx);
    MPI_Reduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        time_stop = MPI_Wtime();
        printf("Parallel time: %3.20f\n", time_stop - time_start);
    }
    return global_max;
}
