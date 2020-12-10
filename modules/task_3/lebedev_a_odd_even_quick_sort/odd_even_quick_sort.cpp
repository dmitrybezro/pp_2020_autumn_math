// Copyright 2020 Lebedev Andrew

#include "../../../modules/task_3/lebedev_a_odd_even_quick_sort/odd_even_quick_sort.h"
#include <mpi.h>
#include <utility>
#include <algorithm>
#include <random>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>


std::vector<int> getRandomVector(int size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> vec(size);
    for (int  i = 0; i < size; i++) {
            vec[i] = gen() % 100;
    }
    return vec;
}


int* odd_even_merge(int* a, int* b, int size1, int size2) {
    int size = size1 + size2;
    int* res = new int[size + 1];
    if (size % 2 == 1) {
        size++;
        res[size - 1] = 1000;
    }
    std::memcpy(&res[0], a, size1 * sizeof(int));
    std::memcpy(&res[0] + size1, b, size2 * sizeof(int));
    for (int p = 1; p < size; p += p)
        for (int k = p; k > 0; k /= 2)
            for (int j = k % p; j + k < size; j += (2 * k))
                for (int i = 0; i < size - j - k; i++)
                    if ((j + i) / (2 * p) == (j + i + k) / (2 * p))
                        if (res[j + i] > res[j + i + k])
                            std::swap(res[j + i], res[j + i + k]);

    return res;
}

int compare(const void *el1, const void *el2) {
  return ( *(static_cast<const int*>(el1)) -
    *(static_cast<const int*>(el2)));
}

int parallel_qsort(int* global_vec, int vec_size) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (vec_size == 0) {
        return 1;
    }
    const int delta = vec_size / size;
    const int rem = vec_size % size;
    int loc_size = delta;
    if (rank < rem) {
        loc_size++;
    }
    if (delta == 0) {
        size = vec_size;
    }
    if (rank == 0) {
        for (int process = 1; process < size; process++) {
            if (process < rem) {
                MPI_Send(global_vec + process * delta + process,
                    delta + 1, MPI_INT, process, 0, MPI_COMM_WORLD);
            } else {
                MPI_Send(global_vec + process * delta + rem,
                    delta, MPI_INT, process, 0, MPI_COMM_WORLD);
            }
        }
    }
    int* local_arr = new int[loc_size];
    if (rank == 0) {
        std::memcpy(local_arr, global_vec, loc_size * sizeof(int));
        qsort(local_arr, loc_size, sizeof(int), compare);
    } else if (rank < size) {
        MPI_Status status;
        MPI_Recv(local_arr, loc_size, MPI_INT,
            0, 0, MPI_COMM_WORLD, &status);
        qsort(local_arr, loc_size, sizeof(int), compare);
    }
    int curr_size = size;
    int level = 1;
    int r_size = 0;
    if (rank < size) {
        while (curr_size > 1) {
            if ( (rank % (2 * level)) == level ) {  // send only
                MPI_Send(&loc_size, 1, MPI_INT, rank - level,
                    0, MPI_COMM_WORLD);
                MPI_Send(local_arr, loc_size, MPI_INT, rank - level,
                    0, MPI_COMM_WORLD);
            }

            if ( !(rank % (2 * level)) && (rank + level < size) ) {
                MPI_Status status;
                MPI_Recv(&r_size, 1, MPI_INT, rank + level, 0,
                    MPI_COMM_WORLD, &status);
                int* r_arr = new int[r_size];
                MPI_Recv(r_arr, r_size, MPI_INT, rank + level, 0,
                    MPI_COMM_WORLD, &status);
                int* tmp = odd_even_merge(local_arr, r_arr, loc_size, r_size);
                loc_size += r_size;
                delete[] r_arr;
                delete[] local_arr;
                local_arr = tmp;
            }

            level *= 2;
            curr_size = curr_size / 2 + curr_size % 2;
        }
    }
    if (rank == 0) {
        std::memcpy(global_vec, local_arr, loc_size * sizeof(int));
    }
    delete[] local_arr;
    return 0;
}
