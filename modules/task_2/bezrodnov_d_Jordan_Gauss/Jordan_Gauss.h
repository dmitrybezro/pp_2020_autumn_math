// Copyright 2020 Bezrodnov Dmitry
#ifndef MODULES_TASK_2_BEZRODNOV_D_JORDAN_GAUSS_JORDAN_GAUSS_H_
#define MODULES_TASK_2_BEZRODNOV_D_JORDAN_GAUSS_JORDAN_GAUSS_H_

#include <vector>

std::vector<double> getRandomMatrix(int size);

std::vector<double> ParallelJordanGauss(const std::vector<double>& MainMatr, int size);

std::vector<double> SequenJordanGauss(const std::vector<double>& MainMatr, int size);

std::vector<double> MultiInverseMatrix(const std::vector<double>& MainMatr, int size);

void FreeMem(double **matr, int n);

void Get_matr(double **matr, int n, double **temp_matr, int indRow, int indCol);

double Det(double **matr, int n);

double** VectorInPointer(const std::vector<double> _matr, int cols);

void Transpon(double **matr, double **tMatr, int n);

std::vector<double> Inverse(double** matr, int n);

std::vector<double> MultiMatrVector(const std::vector<double>& matr, const std::vector<double>& vec);

#endif  //  MODULES_TASK_2_BEZRODNOV_D_JORDAN_GAUSS_JORDAN_GAUSS_H_
