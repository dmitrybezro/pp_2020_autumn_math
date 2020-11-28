// Copyright 2020 Mikhail Vodeneev
#include <mpi.h>
#include <random>
#include <ctime>
#include <vector>
#include <algorithm>
#include "../../../modules/test_tasks/test_mpi/ops_mpi.h"

std::vector<std::vector<double>> getRandomMatrix(int m, int n) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<std::vector<double>> a;
    a.assign(m, std::vector<double>(n));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            a[i][j] = gen() % 10;
    }
    return a;
}

std::vector<double> getRandomMatrixInVector(int size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<double> res(size);
    for (int i = 0; i < size; i++) {
        res[i] = gen() % 100;
    }
    return res;
}

void transpose(std::vector<std::vector<double>>* a) {
    int m = a->size();
    int n = (*a)[0].size();
    std::vector<std::vector<double>> res;
    res.assign(n, std::vector<double>(m));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            res[j][i] = (*a)[i][j];
    }
    *a = res;
}

std::vector<double> matrix_to_vector(std::vector<std::vector<double>> a,
                                                        int m, int n) {
    std::vector<double> res(m * n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            res[i * n + j] = a[i][j];
        }
    }
    return res;
}

std::vector<double> max_el_in_dif_intervals_in_vector(std::vector<double> a,
    int n) {
    double max;
    int delta = a.size() / n;
    std::vector<double> res(delta);
    for (int i = 0; i < delta; i++) {
        max = a[i * n];
        res[i] = max;
        for (int j = 0; j < n; j++) {
            if (a[i * n + j] > max) {
                max = a[i * n + j];
                res[i] = max;
            }
        }
    }
    return res;
}

std::vector<double> getSeqOperations(std::vector<std::vector<double>> a) {
    transpose(&a);
    std::vector<double> res = matrix_to_vector(a, a.size(), a[0].size());
    std::vector<double> globalres = max_el_in_dif_intervals_in_vector(res,
        a[0].size());
    return globalres;
}

std::vector<double> getParallelOperations(std::vector<std::vector<double>> a,
    int m, int n) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int delta = n / size;
    int delta_over = n % size;
    std::vector<double> res;
    if (rank == 0)  {
        transpose(&a);
        res = matrix_to_vector(a, n, m);
        for (int i = 1; i < size; i++)  {
            MPI_Send(&res[0] + i * m * delta,
                m * delta, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    }
    std::vector<double> local_vec(delta*m);
    std::vector<double> local_vec_end(delta_over * m);
    if (rank == 0) {
        local_vec = std::vector<double>(res.begin(), res.begin() + delta*m);
        if (size*delta != n)
            local_vec_end = std::vector<double>(res.end() - delta_over*m,
                res.end());
    } else {
        MPI_Status status;
        MPI_Recv(local_vec.data(), delta*m, MPI_DOUBLE,
            0, 0, MPI_COMM_WORLD, &status);
    }
    std::vector<double> local_res(delta);
    std::vector<double> local_res_end(delta_over);
    local_res = max_el_in_dif_intervals_in_vector(local_vec, m);
    if ((size * delta != n) && rank == 0)
        local_res_end = max_el_in_dif_intervals_in_vector(local_vec_end, m);
    std::vector<double> globalres(n);
    MPI_Gather(local_res.data(), delta, MPI_DOUBLE,
        globalres.data(), delta, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if ((size * delta != n) && rank == 0) {
        for (int i = 0; i < delta_over; i++)
        globalres[n-delta_over+i] = local_res_end[i];
        return globalres;
    }
    return globalres;
}
