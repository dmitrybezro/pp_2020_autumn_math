// Copyright 2018 Nesterov Alexander
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_1/chernyh_d_sum_of_values_ny_matrix_rows/matrix_sum.h"

std::vector<int> getRandomMat(int count_row, int count_str) {
  std::mt19937 gen;
  gen.seed(static_cast<unsigned int>(time(0)));
  std::vector<int> mat(count_row*count_str);
  for (int i = 0; i < count_row*count_str; i++) { mat[i] = gen() % 100; }
  return mat;
}

std::vector<int> getSequentialOperations(std::vector<int> mat, int count_row, int count_str) {
  std::vector<int> reduction_elem(count_str, 0);
  for (int i = 0; i < count_str; i++) {
  for (int j = 0; j < count_row; j++) {reduction_elem[i] += mat[i*count_row + j];}
  }
  return reduction_elem;
}

std::vector<int> getParallelOperations(std::vector<int> global_mat, int count_row, int count_str) {
  int size, size_m, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  size_m = size;
  int delta = count_str / size;
  int outdelta = count_str % size;
  if (size > count_str) {
    size_m = count_str;
    delta = 1;
    outdelta = 0;
  }
  if (rank == 0) {
    for (int proc = 1; proc < size; proc++) {
      if (outdelta != 0 && proc == size_m - 1) {  //  last rank take ostatok
        MPI_Send(&global_mat[0] + proc * delta*count_row, (delta + outdelta)*count_row,
          MPI_INT, proc, 0, MPI_COMM_WORLD);
  } else {
      if (proc >= count_str) {
        MPI_Send(&global_mat[0] + (count_str - 1) * delta*count_row, delta*count_row,
          MPI_INT, proc, 0, MPI_COMM_WORLD);
  } else {
    MPI_Send(&global_mat[0] + proc * delta*count_row, delta*count_row,
      MPI_INT, proc, 0, MPI_COMM_WORLD);
  }
  }
  }
  }
  std::vector<int> local_mat(delta*count_row, 0);
  if (rank == 0) {
    for (int j = 0; j < delta*count_row; j++) { local_mat[j] = global_mat[j]; }
  } else {
    MPI_Status status;
    if (outdelta != 0 && rank == size_m - 1) {
      local_mat.resize((delta + outdelta)*count_row);
      MPI_Recv(&local_mat[0], (delta + outdelta)*count_row, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
  } else {
    MPI_Recv(&local_mat[0], delta*count_row, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
  }
  }
  std::vector<int> global_sum(count_str, 0);
  std::vector<int>local_sum = getSequentialOperations(local_mat, count_row, local_mat.size() / count_row);
  int *recvcount = new int[size];
  int *displs = new int[size];
  for (int i = 0; i < size; i++) {
    if (i >= count_str) {
      recvcount[i] = delta;
      displs[i] = (count_str - 1)*delta;
  } else {
      if (outdelta != 0 && i == size_m - 1) {
        recvcount[i] = delta + outdelta; }
      else
        recvcount[i] = delta;
    displs[i] = i * delta;
  }
  }
  MPI_Gatherv(&local_sum[0], local_sum.size(), MPI_INT, &global_sum[0], recvcount, displs, MPI_INT, 0, MPI_COMM_WORLD);
  delete[]recvcount;
  delete[]displs;
  return global_sum;
}
