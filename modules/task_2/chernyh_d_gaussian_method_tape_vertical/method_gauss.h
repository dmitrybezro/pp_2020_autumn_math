// Copyright 2018 Nesterov Alexander
#ifndef MODULES_TASK_2_CHERNYH_D_GAUSSIAN_METHOD_TAPE_VERTICAL_METHOD_GAUSS_H_
#define MODULES_TASK_2_CHERNYH_D_GAUSSIAN_METHOD_TAPE_VERTICAL_METHOD_GAUSS_H_

#include <vector>
#include <string>

std::vector<double> getRandomMat(int  count_row, int count_str);
std::vector<double> getRandomRes(int count_str);
std::vector<double> Transpose(std::vector<double> A, int count_row, int count_str);
void MatrixPermut(std::vector<double>* T, std::vector<double>* b, int count_row, int count_str);
void MatrixTransform(std::vector<double>* T, std::vector<double>* b, int count_row, int count_str);
double SolutionCheck(std::vector<double> A, std::vector<double> b, std::vector<double> x, int count_row, int count_str);
std::vector<double> getParallelMethod(std::vector<double> global_mat, std::vector<double> b,
  int count_row, int count_str);
std::vector<double> getSequentialMethod(std::vector<double> A, std::vector<double> b, int count_row, int count_str);

#endif  // MODULES_TASK_2_CHERNYH_D_GAUSSIAN_METHOD_TAPE_VERTICAL_METHOD_GAUSS_H_
