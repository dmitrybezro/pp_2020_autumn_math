// Copyright 2020 Mishin Ilya
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_1/mishin_i_matrix_elem_sum/matrix_elem_sum.h"

std::vector<int> getRandomVector(int sz) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> vec(sz);
    for (int i = 0; i < sz; i++) { vec[i] = gen() % 100; }
    return vec;
}

int getSequentialOperations(std::vector<int> vec, std::string ops) {
    const int  sz = vec.size();
    int reduction_elem = 0;
    if (ops == "+") {
        for (int i = 0; i < sz; i++) {
            reduction_elem += vec[i];
        }
    } else if (ops == "-") {
        for (int i = 0; i < sz; i++) {
            reduction_elem -= vec[i];
        }
    } else if (ops == "max") {
        reduction_elem = vec[0];
        for (int i = 1; i < sz; i++) {
            reduction_elem = std::max(reduction_elem, vec[i]);
        }
    }
    return reduction_elem;
}

int getSequentialOperationsMatrix(std::vector<int> vec, int rows, std::string ops) {
    const int  sz = vec.size();
    int reduction_elem = 0;
    if (ops == "+") {
        for (int i = 0; i < sz / rows; i++) {
            for (int j = 0; j < rows; j++) {
                reduction_elem += vec[i * rows + j];
            }
        }
        for (int i = 0; i < sz % rows; i++) {
            reduction_elem += vec[sz - i - 1];
        }
    } else if (ops == "-") {
        for (int i = 0; i < sz / rows; i++) {
            for (int j = 0; j < rows; j++) {
                reduction_elem -= vec[i * rows + j];
            }
        }
        for (int i = 0; i < sz % rows; i++) {
            reduction_elem -= vec[sz - i - 1];
        }
    } else if (ops == "max") {
        reduction_elem = vec[0];
        for (int i = 0; i < sz / rows; i++) {
            for (int j = 0; j < rows; j++) {
                reduction_elem = std::max(reduction_elem, vec[i * rows + j]);
            }
        }
        for (int i = 0; i < sz % rows; i++) {
            reduction_elem = std::max(reduction_elem, vec[sz - i - 1]);
        }
    }
    return reduction_elem;
}


int getParallelOperations(std::vector<int> global_vec,
    int count_size_rows, int count_size_col, std::string ops) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int delta = count_size_rows * count_size_col / size;

    if (rank == 0) {
        for (int proc = 1; proc < size; proc++) {
            MPI_Send(&global_vec[0] + proc * delta, delta,
                MPI_INT, proc, 0, MPI_COMM_WORLD);
        }
    }

    std::vector<int> local_vec(delta);
    if (rank == 0) {
        local_vec = std::vector<int>(global_vec.begin(),
            global_vec.begin() + delta);
    } else {
        MPI_Status status;
        MPI_Recv(&local_vec[0], delta, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }

    int global_sum = 0;
    int local_sum = getSequentialOperationsMatrix(local_vec, count_size_rows, ops);
    MPI_Op op_code;
    if (ops == "+" || ops == "-") { op_code = MPI_SUM; }
    if (ops == "max") { op_code = MPI_MAX; }
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, op_code, 0, MPI_COMM_WORLD);
    return global_sum;
}
