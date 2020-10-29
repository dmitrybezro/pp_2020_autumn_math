// Copyright 2020 Kurnikova Anastasia
#ifndef MODULES_TASK_1_KURNIKOVA_A_INTEGRATION_OPS_MPI_H_
#define MODULES_TASK_1_KURNIKOVA_A_INTEGRATION_OPS_MPI_H_

double f(double x);
int getSequentialOperations(double function(double x),
      double a, double b, int n);
int getParallelOperations(double function(double x),
      double a, double b, int n);

#endif  // MODULES_TASK_1_KURNIKOVA_A_INTEGRATION_OPS_MPI_H_
