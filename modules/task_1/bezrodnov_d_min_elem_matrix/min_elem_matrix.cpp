// Copyright 2020 Bezrodnov Dmitry
#include<mpi.h>
#include <random>
#include <ctime>
#include <vector>
#include <algorithm>
#include "../../../modules/task_1/bezrodnov_d_min_elem_matrix/min_elem_matrix.h"

std::vector<int> getRandomMatrix(int row, int col) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    const int size = row * col;
    std::vector<int> matr(size);
    for (int i = 0; i < row*col; i++) {
        matr[i] = gen() % 1000-i;;
    }
    return matr;
}

int getSequentialOperations(std::vector<int> Matr) {
    int minimum = Matr[0];
    int size = Matr.size();
    for (int i = 0; i < size; i++) {
        minimum = std::min(Matr[i], minimum);
    }
    return minimum;
}


int getParallelOperations(std::vector<int> Matr, int rows, int cols) {
    int RANK;
    int SIZE;

    MPI_Comm_size(MPI_COMM_WORLD, &SIZE);
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);

    int min_loc;
    int min_res;


    const int div = (rows * cols) / SIZE;
    const int mod = (rows * cols) % SIZE;

    if (RANK == 0) {
        for (int i = 1; i < SIZE; i++) {
            if (i <= mod) {
                MPI_Send(&Matr[0] + (div + 1) * i - 1, div + 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            } else {
                MPI_Send(&Matr[0] + div * i + mod, div, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
    }


    std::vector<int> part_vec(div);
    if (RANK == 0) {
        part_vec = std::vector< int>(Matr.begin(),
            Matr.begin() + div);
    } else {
        MPI_Status stat;
        if (RANK <= mod) {
            part_vec.resize(div + size_t(1));
            MPI_Recv(&part_vec[0], div + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
        } else {
            MPI_Recv(&part_vec[0], div, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
        }
    }
    min_loc = getSequentialOperations(part_vec);
    MPI_Reduce(&min_loc, &min_res, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

    return min_res;
}


