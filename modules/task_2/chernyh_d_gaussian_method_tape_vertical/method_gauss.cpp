// Copyright 2018 Nesterov Alexander
#include <mpi.h>
#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <ctime>
#include <algorithm>
#include <utility>
#include "../../../modules/task_2/chernyh_d_gaussian_method_tape_vertical/method_gauss.h"

std::vector<double> getRandomMat(int count_row, int count_str) {
  std::mt19937 gen;
  gen.seed(static_cast<unsigned int>(time(0)));
  std::vector<double> A(count_row*count_str);
  for (int i = 0; i < count_row*count_str; i++) { A[i]= gen() % 9 + 3; }
  return A;
}
std::vector<double> getRandomRes(int count_str) {
  std::mt19937 gen;
  gen.seed(static_cast<unsigned int>(time(0)));
  std::vector<double> b(count_str);
  for (int i = 0; i < count_str; i++) { b[i]= gen() % 9 + 1; }
  return b;
}
std::vector<double> Transpose(std::vector<double> A, int count_row, int count_str) {
  std::vector<double> T(count_row*count_str);
  for (int i = 0; i < count_row; i++)
    for (int j = 0; j < count_str; j++) {
      T[j + i * count_str] = A[i + j * count_row];
    }
  return T;
}

void MatrixPermut(std::vector<double>* T, std::vector<double>* b, int count_row, int count_str) {
  double *pT, *pb;
  pT = T->data();
  pb = b->data();
  double max = pT[0];
  for (int j = 0; j < count_str; j++) {
    max = fabs(pT[j * count_row + j]);
  for (int i = j + 1; i < count_str; i++) {
  if (max < fabs(pT[i * count_row + j])) {
    max = fabs(pT[i * count_row + j]);
    for (int l = 0; l < count_row; l++) {
      std::swap(pT[j * count_row + l], pT[i * count_row + l]);
    }
     std::swap(pb[j], pb[i]); }
  }
  }
  std::vector<double> v1(pT, pT + (count_str*count_row * sizeof(double)) / sizeof(double));
  std::vector<double> v2(pb, pb + (count_str * sizeof(double)) / sizeof(double));
  T = &v1;
  b = &v2;
}

void MatrixTransform(std::vector<double>* T, std::vector<double>* b, int count_row, int count_str) {
  double *pT, *pb;
  pT = T->data();
  pb = b->data();
  double kf = 0;
  for (int i = 0; i < count_str; i++) {
    for (int j = i + 1; j < count_str; j++) {
      if (pT[i*count_row + i] == 0) {
        while (pT[i*count_row + i] == 0 && i < count_str - 1) {
          for (int l = 0; l < count_row; l++) {
            std::swap(pT[i * count_row + l], pT[(i+1) * count_row + l]);
          }
          std::swap(pb[i], pb[i+1]);
        }
      } else {
        kf = pT[j*count_row + i] / pT[i*count_row + i];
        if (kf != 0) {
          for (int k = i; k < count_row; k++) {
            pT[j*count_row + k] -= pT[i*count_row + k] * kf;
          }
          pb[j] -= pb[i] * kf;
        }
      }
    }
  }
  std::vector<double> v1(pT, pT + (count_str*count_row*sizeof(double)) / sizeof(double));
  std::vector<double> v2(pb, pb + (count_str*sizeof(double)) / sizeof(double));
  T = &v1;
  b = &v2;
}
double SolutionCheck(std::vector<double> A, std::vector<double> b,
    std::vector<double> x, int count_row, int count_str) {
  std::vector<double> error(count_str, 0);
  double max_error = 0;
  for (int i = 0; i < count_str; i++) {
    for (int j = 0; j < count_row; j++) {
      error[i] += A[i*count_row + j] * x[j];
    }
    error[i] -= b[i];
  }
  max_error = error[0];
  for (int i = 0; i < count_str; i++) {
  if (max_error < error[i]) max_error = error[i];
  }
  return max_error;
}


std::vector<double> getSequentialMethod(std::vector<double> A, std::vector<double> b, int count_row, int count_str) {
  std::vector<double> x(count_row, 0);
  double sum = 0;
  if (count_str != count_row || count_str <= 0) {
    throw std::runtime_error("Wrong size");
  }
  MatrixPermut(&A, &b, count_row, count_str);
  // calculate coeff
  MatrixTransform(&A, &b, count_row, count_str);
  MatrixPermut(&A, &b, count_row, count_str);
  // calculate x
  for (int i = count_row - 1; i >= 0; i--) {
    sum = 0;
    if (b[i] != 0 && A[i*count_row + i] == 0) {
      throw std::runtime_error("Inconistent system");
    } else {
        for (int j = i; j < count_str; j++) {
          sum += x[j] * A[i*count_row + j];
          x[i] = (b[i] - sum) / A[i*count_row + i];
        }
    }
  }
  return x;
}

std::vector<double> getParallelMethod(std::vector<double> global_mat, std::vector<double> b,
    int count_row, int count_str) {
  int size, rank, size_m, code = 0;
  std::vector<int> Iter(count_str, -1);
  std::vector<int> LeadPos(count_str, -1);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  size_m = size;
  int delta = count_row / size;
  int outdelta = count_row % size;
  if (count_str != count_row || count_str <= 0) {
    code = -1;
    MPI_Bcast(&code, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  if (code != 0) {
    throw std::runtime_error("Wrong size");
  }
  std::vector<double> global_mat_t(count_row*count_str, 0);
  if (size > count_row) {
    size_m = count_row;
    delta = 1;
    outdelta = 0;
  }
  if (rank == 0) {
    global_mat_t = Transpose(global_mat, count_row, count_str);
    for (int proc = 1; proc < size; proc++) {
      if (outdelta != 0 && proc == size_m - 1) {
        MPI_Send(&global_mat_t[0] + proc * delta*count_str, (delta + outdelta)*count_str,
          MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
      } else {
          if (proc >= count_row) {
            MPI_Send(&global_mat_t[0] + (count_row - 1) * delta*count_str, delta*count_str,
              MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
          } else {
              MPI_Send(&global_mat_t[0] + proc * delta*count_str, delta*count_str,
                MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
          }
      }
    }
  }
  MPI_Bcast(&b[0], count_str, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  std::vector<double> local_mat(delta*count_str, 0);
  if (rank == 0) {
    for (int j = 0; j < delta*count_str; j++) { local_mat[j] = global_mat_t[j]; }
  } else {
      MPI_Status status;
      if (outdelta != 0 && rank == size_m - 1) {
        local_mat.resize((delta + outdelta)*count_str);
        MPI_Recv(&local_mat[0], (delta + outdelta)*count_str, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
      } else {
          MPI_Recv(&local_mat[0], delta*count_str, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
      }
  }
  for (int proc = 0; proc < size_m; proc++) {
    double max = 0;
    int pos_max = -1;
    if (proc == rank) {
      int k = rank * delta;
      for (int i = 0; i < static_cast<int>(local_mat.size() / count_str); i++) {
        max = 0;
        for (int j = 0; j < count_str; j++) {
          if (Iter[j] == -1 && max < fabs(local_mat[i*count_str + j])) {
            max = fabs(local_mat[i*count_str + j]);
            pos_max = j;
          }
        }
        if (Iter[pos_max] == -1 && max == fabs(local_mat[i*count_str + pos_max])) {
          Iter[pos_max] = k + i;
          LeadPos[Iter[pos_max]] = pos_max;
        }
      }
    }
    MPI_Bcast(&Iter[0], count_str, MPI_INT, proc, MPI_COMM_WORLD);
    MPI_Bcast(&LeadPos[0], count_str, MPI_INT, proc, MPI_COMM_WORLD);
  }
  std::vector<double> kf(count_str, 0);
  int flag = 0, lp = -1, f_swap = 0;
  int delta_flag = 0;
  for (int proc = 0; proc < size_m; proc++) {
    if (proc == rank) {
      f_swap = 0;
      delta = local_mat.size() / count_str;
      delta_flag = (rank + 1)*delta;
      if (outdelta != 0 && rank == size_m - 1) { delta_flag = count_str; }
      for (int proc_t = 0; proc_t < size_m; proc_t++) {
        if (proc_t != proc)
          MPI_Send(&delta, 1, MPI_INT, proc_t, 0, MPI_COMM_WORLD);
      }
      if (flag < delta_flag) {
        for (int i = 0; i < static_cast<int>(local_mat.size() / count_str); i++) {
          lp = LeadPos[flag];
          if (local_mat[i*count_str + lp] == 0) {
            f_swap = 1;
            while (local_mat[i*count_str + lp] == 0 && flag < count_str - 1) {
              std::swap(LeadPos[flag], LeadPos[flag + 1]);
              lp = LeadPos[flag];
            }
            for (int proc_t = 0; proc_t < size_m; proc_t++) {
              if (proc_t != proc)
                MPI_Send(&LeadPos[0], count_str, MPI_DOUBLE, proc_t, 3, MPI_COMM_WORLD);
            }
          }
          for (int p = 0; p < count_str; p++) kf[p] = 0;
          for (int j = 0; j < count_str; j++) {
            if (j != lp && local_mat[i*count_str + lp] != 0) {
              kf[j] = local_mat[i*count_str + j] / local_mat[i*count_str + lp];
              if (kf[j] != 0) {
                for (int s = i; s < static_cast<int>(local_mat.size() / count_str); s++) {
                  local_mat[s*count_str + j] -= local_mat[s*count_str + lp] * kf[j];
                }
                b[j] -= b[lp] * kf[j];
              }
            }
          }
          for (int proc_t = 0; proc_t < size_m; proc_t++) {
            if (proc_t != proc) {
              MPI_Send(&lp, 1, MPI_INT, proc_t, 1, MPI_COMM_WORLD);
              MPI_Send(&kf[0], count_str, MPI_DOUBLE, proc_t, 2, MPI_COMM_WORLD);
            }
          }
          flag++;
        }
      }
    }
    MPI_Bcast(&f_swap, 1, MPI_INT, proc, MPI_COMM_WORLD);
    if (proc != rank && rank < size_m) {
      int num = 0, f = 0;
      MPI_Status status0, status1, status2;
      MPI_Recv(&num, 1, MPI_INT, proc, 0, MPI_COMM_WORLD, &status0);
      if (f_swap) {
        MPI_Status status;
        MPI_Recv(&LeadPos[0], count_str, MPI_DOUBLE, proc, 3, MPI_COMM_WORLD, &status);
      }
      while (f < num) {
        MPI_Recv(&lp, 1, MPI_INT, proc, 1, MPI_COMM_WORLD, &status1);
        MPI_Recv(&kf[0], count_str, MPI_DOUBLE, proc, 2, MPI_COMM_WORLD, &status2);
        for (int i = 0; i < static_cast<int>(local_mat.size() / count_str); i++) {
          for (int j = 0; j < count_str; j++) {
            if (j != lp) {
              if (kf[j] != 0) {
                local_mat[i*count_str + j] -= local_mat[i*count_str + lp] * kf[j];
              }
            }
          }
        }
        f++;
      }
    }
    MPI_Bcast(&flag, 1, MPI_INT, proc, MPI_COMM_WORLD);
    MPI_Bcast(&b[0], count_str, MPI_DOUBLE, proc, MPI_COMM_WORLD);
  }
  if (rank == size_m - 1 && size_m < size) {
    for (int p = size_m; p < size; p++) {
      MPI_Send(&local_mat[0], delta*count_str, MPI_DOUBLE, p, 4, MPI_COMM_WORLD);
    }
  }
  if (rank >= size_m) {
    MPI_Status status;
    MPI_Recv(&local_mat[0], delta*count_str, MPI_DOUBLE, size_m - 1, 4, MPI_COMM_WORLD, &status);
  }
  std::vector<double> rez(local_mat.size() / count_str, 0);
  int flag_rez = 0, lp_rez = -1;
  for (int proc = 0; proc < size; proc++) {
    if (proc == rank && proc < size_m) {
      delta_flag = (rank + 1)*delta;
      if (outdelta != 0 && rank == size_m - 1) {
        delta_flag = count_str;
      }
      if (flag_rez < delta_flag) {
        for (int i = 0; i < static_cast<int>(rez.size()); i++) {
          lp_rez = LeadPos[flag_rez];
          if (b[lp_rez] != 0 && local_mat[i*count_str + lp_rez] == 0) {
            code = -1;
            break;
          } else {
              rez[i] = b[lp_rez] / local_mat[i*count_str + lp_rez];
          }
          flag_rez++;
        }
      }
    }
    if (rank >= size_m && code != -1) {
      lp_rez = LeadPos[count_str - 1];
      rez[0] = b[lp_rez] / local_mat[lp_rez];
    }
    MPI_Bcast(&flag_rez, 1, MPI_INT, proc, MPI_COMM_WORLD);
    MPI_Bcast(&code, 1, MPI_INT, proc, MPI_COMM_WORLD);
  }
  if (code == -1) {
    throw std::runtime_error("Inconistent system");
  }
  std::vector<double> global_rez(count_row, 0);
  std::vector<int> recvcount(size, 0);
  std::vector<int> displs(size, 0);
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
  MPI_Gatherv(&rez[0], rez.size(), MPI_DOUBLE, &global_rez[0], &recvcount[0],
               &displs[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return global_rez;
}
