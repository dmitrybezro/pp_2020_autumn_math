// Copyright 2020 Vodeneev Mikhail
#include <vector>
#include <string>

std::vector<std::vector<double>> getRandomMatrix(int m, int n);

std::vector<double> getRandomMatrixInVector(int size);

void transpose(std::vector<std::vector<double>>* a);

std::vector<double> matrix_to_vector(std::vector<std::vector<double>> a,
     int m, int n);

std::vector<double> max_el_in_dif_intervals_in_vector(std::vector<double> a,
     int n);

std::vector<double> getSeqOperations(std::vector<std::vector<double>> a);

std::vector<double> getParallelOperations(std::vector<std::vector<double>> a,
     int m, int n);
