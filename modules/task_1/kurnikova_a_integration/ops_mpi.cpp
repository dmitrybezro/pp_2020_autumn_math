// Copyright 2020 Kurnikova Anastasia
#include <mpi.h>
#include "../../../modules/task_1/kurnikova_a_integration/ops_mpi.h"

double f(double x) {
    return (4 / (1 + x * x));
}

int getSequentialOperations(double function(double x),
      double a, double b, int n) {
    const double h = (b - a) / static_cast<double>(n);
    double answer = 0;
    for (int i = 0; i < n; i++)
        answer = answer + function(a + h * (i + 0.5));
    answer = h * answer;
    return answer;
}

int getParallelOperations(double function(double x),
      double a, double b, int n) {
    int procnum, procrank;
    const double h = (b - a) / static_cast<double>(n);
    double answer = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &procnum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    double procans = 0;
    for (int i = procrank; i < n; i += procnum)
        procans = procans + function(a + h * (i + 0.5));
    MPI_Reduce(&procans, &answer, 1, MPI_DOUBLE,
      MPI_SUM, 0, MPI_COMM_WORLD);
    answer = h * answer;
    return answer;
}
