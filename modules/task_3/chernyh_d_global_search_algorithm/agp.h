// Copyright 2018 Nesterov Alexander
#ifndef MODULES_TASK_3_CHERNYH_D_GLOBAL_SEARCH_ALGORITHM_AGP_H_
#define MODULES_TASK_3_CHERNYH_D_GLOBAL_SEARCH_ALGORITHM_AGP_H_

#include <vector>
#include <string>

std::vector<int> getRandomVector(int  sz);
void CoefIn(std::vector<int> _coef);
void FuncInit(void);
double Func(int index_funk, double x);
struct Trial {
  double x, z;
  Trial& operator = (double val) {
    x = val;
    z = val;
    return *this;
  }
  Trial& operator = (const Trial &tr) {
    x = tr.x;
    z = tr.z;
    return *this;
  }
  bool operator == (const Trial &tr) {
    double epsilon = 0.00001;
    if (fabs(x - tr.x) < epsilon && fabs(z - tr.z) < epsilon)
     return true;
    else
     return false;
  }
};
Trial getSequential(std::vector<Trial> trials, int index_func, double eps, double r);
Trial getParallelOperations(std::vector<Trial> trials, int index_func, double eps, double r);

#endif  // MODULES_TASK_3_CHERNYH_D_GLOBAL_SEARCH_ALGORITHM_AGP_H_
