// Copyright 2020 Kolesnikov Gleb
#ifndef MODULES_TASK_3_KOLESNIKOV_G_TRAPEZOID_RULE_OPS_MPI_H_
#define MODULES_TASK_3_KOLESNIKOV_G_TRAPEZOID_RULE_OPS_MPI_H_

#include <mpi.h>
#include <vector>

double trapezoidParallelRule(double(*func)(std::vector<double>), std::vector <double> start_coords,
    std::vector <double> end_coords, unsigned int n);
double SequentialIntegr_smart(double(*func)(std::vector<double>), std::vector <double> x,
    std::vector <double> y, const int n);

#endif  // MODULES_TASK_3_KOLESNIKOV_G_TRAPEZOID_RULE_OPS_MPI_H_
