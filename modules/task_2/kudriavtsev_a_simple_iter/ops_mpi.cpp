// Copyright 2020 Kudriavtsev Alexander
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <iostream>
#include "./ops_mpi.h"

std::vector<double> getRandomVector(int n) {
    std::vector<double> vec(n);
    double lower_bound = -100.0;  // Why?
    double upper_bound = 100.0;  // Why?
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    std::default_random_engine re;
    for (int i = 0; i < n; ++i) {
        vec[i] = unif(re);
    }
    return vec;
}

std::vector<double> getRandomMatrix(int n) {
    std::vector<double> vec(n * n);
    double lower_bound = -0.99;
    double upper_bound = 0.99;
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    std::default_random_engine re;
    for (int i = 0; i < n; ++i) {  // Could be done in one "for", but who cares?
        for (int j = 0; j < n; ++j) {
            vec[i*n+j] = unif(re);
        }
        vec[i * n + i] = n + 2;
    }
    return vec;
}

std::vector<double> parallelIterMethod(std::vector<double> A, std::vector<double> b, int n, double eps, int NMax) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int delta = n / size;
    const int rem = n % size;  // In the first lab I made a mistake. REMAINDER not REMEDIAL. heh.
    double t_b, t_e;
    std::vector<double> x(n);
    if (n < size) {
        return (rank == 0) ?  sequentialIterMethod(A, b, n, eps, NMax): std::vector<double>();
    }
    if (rank == 0) {
        t_b = MPI_Wtime();
    }
    // -- Things for Scatterv and Allgatherv --
    int* counts = new int[size * 4];  // ONE new/delete is better than FOUR. HEH.
    int* displs_b = counts, *displs_A = counts + size;
    int *sendcounts_b = counts + size + size, *sendcounts_A = counts + size + size + size;
    displs_A[0] = displs_b[0] = 0;
    sendcounts_b[0] = delta + static_cast<int>(0 < rem);
    sendcounts_A[0] = (delta + static_cast<int>(0 < rem)) * n;
    for (int i = 1; i < size; ++i) {
        displs_b[i] = displs_b[i-1] + delta + static_cast<int>(i <= rem);
        displs_A[i] = displs_A[i - 1] + (delta + static_cast<int>(i <= rem)) * n;
        sendcounts_b[i] = delta + static_cast<int>(i < rem);
        sendcounts_A[i] = (delta + static_cast<int>(i < rem)) * n;
    }
    // -- Data for every procces --
    int proc_size = delta + static_cast<int>(rank < rem);
    double* A_proc = new double[proc_size * n];
    double* b_proc = new double[proc_size];
    double* x_proc = new double[proc_size];
    double* tmp = new double[n];
    MPI_Scatterv(A.data(), sendcounts_A, displs_A, MPI_DOUBLE, A_proc, proc_size * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {  // Preparation. Needs to be done in root process.
        for (int i = 0; i < n; ++i) {
            x[i] = tmp[i] = b[i] = b[i] / A[i * n + i];
        }
    }
    MPI_Scatterv(b.data(), sendcounts_b, displs_b, MPI_DOUBLE, b_proc, proc_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(tmp, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // -- Actual Algorithm --
    for (int i = 0; i < proc_size; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j != i)
                A_proc[i * n + j] = -(A_proc[i * n + j] / A_proc[i * n + i]);
        }
        A_proc[i * n + i] = 0.0;
    }
    for (int i = 0; i < n; ++i) {
        x[i] = tmp[i];
    }
    double max_diff = 0.0;
    int iter_count = 0;
    do {
        ++iter_count;
        for (int i = 0; i < proc_size; ++i) {
            x_proc[i] = b_proc[i];
            for (int j = 0; j < n; ++j) {
                x_proc[i] = x_proc[i] + (A_proc[i * n + j] * tmp[j]);
            }
        }
        MPI_Allgatherv(x_proc, proc_size, MPI_DOUBLE, x.data(), sendcounts_b, displs_b, MPI_DOUBLE, MPI_COMM_WORLD);
        max_diff = 0.0;
        double local_md = 0.0;
        for (int i = 0; i < proc_size; ++i) {  // Could be done in one process or in all processes
            if (fabs(x_proc[i] - tmp[displs_b[rank] + i]) > local_md) {
                local_md = fabs(x_proc[i] - tmp[displs_b[rank] + i]);
            }
        }
        MPI_Allreduce(&local_md, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        for (int i = 0; i < n; ++i) {
            tmp[i] = x[i];
        }
    } while (max_diff >= eps && iter_count < NMax);
    // -- Release RECOLLEC..HM! MEMORY --
    delete[] A_proc;  // IDK Why everyone always forgets about delete.
    delete[] b_proc;
    delete[] x_proc;
    delete[] counts;
    delete[] tmp;
    if (rank == 0) {
        t_e = MPI_Wtime();
        std::cout << "Parallel time: " << t_e - t_b << std::endl;
    }
    return x;
}

std::vector<double> sequentialIterMethod(std::vector<double> A, std::vector<double> b, int n, double eps, int NMax) {
    double t_b, t_e;  // Time of the Beginning and of the End...
    std::vector<double> x(n);
    double* tmp_arr = new double[n];
    t_b = MPI_Wtime();
    for (int i = 0; i < n; ++i) {  // Preparations
        tmp_arr[i] = b[i] / A[i * n + i];
        b[i] = b[i] / A[i * n + i];
        for (int j = 0; j < n; ++j) {
            if (j != i)
                A[i * n + j] = -(A[i * n + j] / A[i * n + i]);
        }
        A[i * n + i] = 0.0;
    }
    int iter_count = 0;
    double maxdiff = 0.0;
    do {  // rare case of using do ... while(...).
        ++iter_count;
        for (int i = 0; i < n; ++i) {  // Preparations
            x[i] = b[i];
            for (int j = 0; j < n; ++j) {
                x[i] += A[i * n + j] * tmp_arr[j];
            }
        }
        maxdiff = 0.0;
        for (int i = 0; i < n; ++i) {
            if (fabs(x[i] - tmp_arr[i]) > maxdiff) {
                maxdiff = fabs(x[i] - tmp_arr[i]);
            }
            tmp_arr[i] = x[i];
        }
    } while (maxdiff >= eps && iter_count < NMax);
    delete[] tmp_arr;
    t_e = MPI_Wtime();
    std::cout << "Sequential Time :" << t_e - t_b << std::endl;
    return x;
}

double discrepancyNorm(const std::vector<double>& x, const std::vector<double>& A, const std::vector<double>& b) {
    const int n = x.size();
    double dis = 0.0;
    for (int i = 0; i < n; ++i) {
        double tmp = 0.0;
        for (int j = 0; j < n; ++j) {
            tmp += A[i * n + j] * x[j];
        }
        if (fabs(tmp-b[i]) > dis) {
            dis = fabs(tmp-b[i]);
        }
    }
    return dis;
}
