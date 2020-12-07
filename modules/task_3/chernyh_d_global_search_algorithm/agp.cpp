// Copyright 2018 Nesterov Alexander
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include<iostream>
#include "../../../modules/task_3/chernyh_d_global_search_algorithm/agp.h"

std::vector<int> coef(5);

void CoefIn(std::vector<int> _coef) {
  coef = _coef;
}

void FuncInit(void) {
  std::mt19937 gen;
  gen.seed(static_cast<unsigned int>(time(0)));
  for (int i = 0; i < coef.size(); i++) coef[i] = gen() % 10 + 1;
}

double Func(int index_funk, double x) {
  switch (index_funk) {
    case 0:
      return fabs(x - 1);
      break;
    case 1:
      return (x*x + 2);
      break;
    case 2:
      return (coef[0] * sin(coef[1] * x) + coef[2] * cos(coef[3] * x) + coef[4]);
      break;
    default:
      throw std::runtime_error("Wrong index of funktion");
  }
}

Trial getSequential(std::vector<Trial> trials, int index_func, double eps, double r) {
  double M, Rmax, Rpos;
  Trial first = trials[0];
  Trial second = trials[1];
  Trial current = trials[0];
  double curr_eps = second.x - first.x;
  std::vector<Trial>::iterator it = trials.begin();
  while (curr_eps > eps) {
  Rpos = 1;
  M = fabs((trials[1].z - trials[0].z) / (trials[1].x - trials[0].x));
  for (int i = 2; i < static_cast<int>(trials.size()); i++) {
    double max;
    max = fabs((trials[i].z - trials[i - 1].z) / (trials[i].x - trials[i - 1].x));
    if (max > M)
      M = max;
  }
  if (M == 0) M = 1;
  else
    M = r * M;
  Rmax = M * (trials[1].x - trials[0].x) + (pow((trials[1].z - trials[0].z), 2)
      / (M * (trials[1].x - trials[0].x)))- 2 * (trials[1].z + trials[0].z);
  for (int i = 2; i <static_cast<int>(trials.size()); i++) {
    double k = M * (trials[i].x - trials[i - 1].x);
    double R = k + (pow((trials[i].z - trials[i - 1].z), 2) / k)
        - 2 * (trials.at(i).z + trials.at(i - 1).z);
    if (R > Rmax) {
      Rmax = R;
      Rpos = i;
    }
  }
  curr_eps = trials.at(Rpos).x - trials.at(Rpos - 1).x;
  std::vector<Trial>::iterator it2 = trials.begin();
  for (it = trials.begin(); it - trials.begin() <= Rpos; it++) it2 = it;
  current.x = (trials.at(Rpos).x + trials.at(Rpos - 1).x) / 2
      - (trials.at(Rpos).z - trials.at(Rpos - 1).z) / (2 * M);
  current.z = Func(index_func, current.x);
  trials.insert(it2, current);
  }
  return current;
}

struct Charact {
  double r;
  int pos;
};

Trial getParallelOperations(std::vector<Trial> trials, int index_func, double eps, double r) {
  int size, size_m, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Bcast(&trials[0], 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  Trial first = trials[0];
  Trial second = trials[1];
  Trial current = trials[0];
  double curr_eps = second.x - first.x;
  std::vector<Trial>::iterator it = trials.begin();
  int itr = 0;
  int best_i = 0;
  double best_z = first.z;
  double Mmax = 0;
  while (curr_eps > eps) {
    size_m = size;
    int trials_size = static_cast<int> (trials.size());
    int delta = (trials_size - 1) / size;
    int outdelta = (trials_size - 1) % size;
    if (size > trials_size - 1) {
      delta = 1;
      outdelta = 0;
      size_m = trials_size - 1;
    }
    int it_start = rank * delta;
    int it_finish = (rank + 1)*delta;
    if (outdelta != 0 && rank == size - 1) it_finish = (rank + 1)*delta + outdelta;
    if (rank >= size_m) {
      it_start = (size_m - 1)* delta;
      it_finish = (size_m)*delta;
    }
    Charact ch_interval;
    ch_interval.pos = it_start+1;
    double M;
    M = fabs((trials[it_start + 1].z - trials[it_start].z) / (trials[it_start + 1].x - trials[it_start].x));
    for (int i = it_start + 2; i <= it_finish; i++) {
      double max;
      max = fabs((trials[i].z - trials[i - 1].z) / (trials[i].x - trials[i - 1].x));
      if (max > M) M = max;
    }
    MPI_Allreduce(&M, &Mmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (Mmax == 0)
      Mmax = 1;
    else
      Mmax = r * Mmax;
    Charact Rmax;
    ch_interval.r = Mmax * (trials[it_start + 1].x - trials[it_start].x) + (pow((trials[it_start + 1].z
       - trials[it_start].z), 2) / (Mmax * (trials[it_start + 1].x - trials[it_start].x)))
       - 2 * (trials[it_start + 1].z + trials[it_start].z);
    for (int i = it_start + 2; i <= it_finish; i++) {
      double k, R;
      k = Mmax * (trials[i].x - trials[i - 1].x);
      R = k + (pow((trials[i].z - trials[i - 1].z), 2) / k) - 2 * (trials[i].z + trials[i - 1].z);
      if (R > ch_interval.r) {
        ch_interval.r = R;
        ch_interval.pos = i;
      }
    }
    MPI_Allreduce(&ch_interval, &Rmax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    curr_eps = trials[Rmax.pos].x - trials[Rmax.pos - 1].x;
    current.x = (trials[Rmax.pos].x + trials[Rmax.pos - 1].x) / 2
              - (trials[Rmax.pos].z - trials[Rmax.pos - 1].z) / (2 * Mmax);
    current.z = Func(index_func, current.x);
    std::vector<Trial>::iterator it2 = trials.begin();
    for (it = trials.begin(); it - trials.begin() <= Rmax.pos; it++) it2 = it;
    trials.insert(it2, current);
  }
  return current;
}
