// Copyright 2020 Kolesnikov Gleb

#include "./ops_mpi.h"
#include <mpi.h>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>


double SequentialIntegr_smart(double(*func)(std::vector<double>), std::vector <double> x,
    std::vector <double> y, const int n) {

    int dimension = 3;
    int borderCounter;
    std::vector<double> h(dimension);

    for (int i = 0; i < dimension; ++i) {
        h[i] = (y[i] - x[i]) / static_cast<double>(n);
    }
    std::vector <double> segments1(dimension);
    double result = 0.0;
    int m = n + 1;
    int num = std::pow(n + 1, dimension);
    for (int i = 0; i < num; ++i) {
        borderCounter = 0;
        segments1[0] = x[0] + h[0] * (i % m);
        segments1[1] = x[1] + h[1] * ((i % (m*m)) / m);
        segments1[2] = x[2] + h[2] * (i / (m*m));
        for (int j = 0; j < dimension; j++) {
            if (segments1[j] == x[j] || segments1[j] == y[j])
                borderCounter++;
        }
        result += func(segments1) / std::pow(2, borderCounter);
    }
    for (int i = 0; i < dimension; ++i) {
        result *= h[i];
    }
    return result;
}

double trapezoidParallelRule(double(*func)(std::vector<double>), std::vector <double> start_coords,
    std::vector <double> end_coords, unsigned int n) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int m = n + 1;
    int dimension = 3;
    int num = std::pow(m, dimension);
    int locstart, locend;
    int delta = num / size;
    int rem = num % size;
    double locres = 0.0;
    double res = 0.0;
    int borderCounter = 0;
    std::vector<double> h(dimension);
    for (int i = 0; i < dimension; ++i) {
        h[i] = (end_coords[i] - start_coords[i]) / static_cast<double>(n);
    }
    if (rank < rem) {
        locstart = rank*(delta + 1);
        locend = rank * (delta + 1) + delta + 1;
    }  else {
        locstart = rank * delta + rem;
        locend = rank * delta + delta + rem;
    }
    std::vector <double> point(dimension);
    for (int i = locstart; i < locend; ++i) {
        borderCounter = 0;
        point[0] = start_coords[0] + h[0] * (i % m);
        point[1] = start_coords[1] + h[1] * ((i % (m*m)) / m);
        point[2] = start_coords[2] + h[2] * (i / (m*m));
        for (int j = 0; j < dimension; j++) {
            if (point[j] == start_coords[j] || point[j] == end_coords[j])
                borderCounter++;
        }
        locres += func(point) / std::pow(2, borderCounter);
    }
    MPI_Reduce(&locres, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    res *= h[0] * h[1] * h[2];
    return res;
}
