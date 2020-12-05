// Copyright 2020 Kurnikova Anastasia
#include <mpi.h>
#include <vector>
#include "../../../modules/task_2/kurnikova_a_jacobi/ops_mpi.h"

int correctness(std::vector<double> a, int n) {
    double abssum = 0;
    for (int i = 0; i < n; i++) {
        if (a[i * n + i] == 0)
            return 0;
        for (int j = 0; j < n; j++)
            if (i != j)
                abssum = abssum + fabs(a[i * n + j] / a[i * n + i]);
        if (abssum > 1)
            return 0;
    }
    return 1;
}

int sequential(std::vector<double> a, std::vector<double> b,
                                   int n, const double e) {
    int it = 0;
    double md = 0;
    std::vector<double> x0(n);
    std::vector<double> tmp(n);
    do {
        for (int i = 0; i < n; i++) {
            tmp[i] = b[i];
            for (int j = 0; j < n; j++)
                if (i != j)
                    tmp[i] = tmp[i] - a[i * n + j] * x0[j];
            tmp[i] = tmp[i] / a[i * n + i];
        }
        for (int k = 0; k < n; k++) {
            x0[k] = b[k];
            if (fabs(x0[k] - tmp[k]) > md)
                md = fabs(x0[k] - tmp[k]);
        }
        it++;
    } while ((md >= e) && (it < 100));
    if (it == 99)
        return 0;
    else
        return 1;
}

int parallel(std::vector<double> a, std::vector<double> b,
                                 int n, const double e) {
    int procnum, procrank, oneprocsize, it = 0;
    double md = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &procnum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    oneprocsize = n / procnum;
    std::vector<double> xproc(n);
    std::vector<double> tmp(n);
    std::vector<double> x0(n);
    std::vector<double> aproc(oneprocsize * n);
    std::vector<double> bproc(oneprocsize);
    MPI_Scatter(a.data(), oneprocsize * n, MPI_DOUBLE,
                aproc.data(), oneprocsize * n, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
    MPI_Scatter(b.data(), oneprocsize, MPI_DOUBLE,
                bproc.data(), oneprocsize, MPI_DOUBLE,
                               0, MPI_COMM_WORLD);
    for (int i = 0; i < oneprocsize; i++)
        xproc[i] = bproc[i];
    do {
        for (int i = 0; i < n; i++)
            tmp[i] = b[i];
        for (int i = 0; i < oneprocsize; i++) {
            xproc[i] = bproc[i];
            for (int j = 0; j < n; j++)
                if ((procrank * oneprocsize + i) != j)
                    xproc[i] = xproc[i] - x0[j] * aproc[i * n + j];
            xproc[i] = xproc[i] / aproc[i * n + procrank * oneprocsize + i];
        }
        MPI_Allgather(xproc.data(), oneprocsize, MPI_DOUBLE,
                      tmp.data(), oneprocsize, MPI_DOUBLE,
                                      MPI_COMM_WORLD);
        for (int k = 0; k < n; k++) {
            x0[k] = b[k];
            if (fabs(x0[k] - tmp[k]) > md)
                md = fabs(x0[k] - tmp[k]);
        }
        it++;
    } while ((md >= e) && (it < 100));
    if (it == 99)
        return 0;
    else
        return 1;
}
