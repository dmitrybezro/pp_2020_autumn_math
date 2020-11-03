  // Copyright 2020 Tsvetkov Maxim
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include <iostream>
#include "../../../modules/task_1/tsvetkov_m_vector_disorderly_order/vector_disorderly_order.h"

std::vector<int> getRandomVector(int sz) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> vec(sz);
    for (int i = 0; i < sz; i++) { vec[i] = gen() % 100; }
    return vec;
}

int get_sequential_operations(std::vector<int> vec) {
    const int  sz = vec.size();
    int disorder = 0;
    for (int i = 0; i < sz - 1; i++) {
        if (vec[i] > vec[i + 1]) {
            ++disorder;
        }
    }

    return disorder;
}

int get_parallel_operations(std::vector<int> global_vec, int vec_size) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int delta = vec_size / size;
    const int k = vec_size % size;

    if (rank == 0) {
        for (int proc = 1; proc < size; proc++) {
            if (proc < k) {
                MPI_Send(&global_vec[0] + proc * delta + (proc - 1) + 1, delta + 1,
                    MPI_INT, proc, 0, MPI_COMM_WORLD);
            } else {
                MPI_Send(&global_vec[0] + proc * delta + k, delta,
                    MPI_INT, proc, 0, MPI_COMM_WORLD);
            }
        }
    }

    std::vector<int> local_vec(delta);
    if (rank == 0) {
        local_vec = std::vector<int>(global_vec.begin(),
            global_vec.begin() + delta + ((k == 0) ? 0 : 1));
    } else {
        MPI_Status status;
        if (rank < k) {
            local_vec.resize(delta + 1);
            MPI_Recv(&local_vec[0], delta + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        } else {
            MPI_Recv(&local_vec[0], delta, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        }
    }


    int global_sum = 0;
    int local_sum = get_sequential_operations(local_vec);

    if (rank == 0) {
        for (size_t i = 1; i < size; i++) {
            if (i < k) {
                if (global_vec[i * (delta + 1) - 1] > global_vec[i * (delta + 1)]) {
                    ++local_sum;
                }
            } else {
                if (global_vec[i * (delta)+k - 1] > global_vec[i * (delta)+k]) {
                    ++local_sum;
                }
            }
        }
    }

    MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_sum;
}

