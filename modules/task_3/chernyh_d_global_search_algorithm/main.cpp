// Copyright 2018 Nesterov Alexander
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include "./agp.h"


TEST(Parallel_Operations_MPI, TestAGP_1) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<Trial> trials(2);
  int in_func = 0;
  double x_0 = 0, x_n = 3, x_true = 1;
  double eps = 0.0001, r = 2.5;
  if (rank == 0) {
    Trial f, s;
    f.x = x_0;
    f.z = Func(in_func, x_0);
    trials[0] = f;
    s.x = x_n;
    s.z = Func(in_func, x_n);
    trials[1] = s;
  }
  double time1 = MPI_Wtime();
  Trial rez1 = getParallelOperations(trials, in_func, eps, r);
  double time2 = MPI_Wtime();
  EXPECT_NEAR(x_true, rez1.x, eps);
  if (rank == 0) {
    double time3 = MPI_Wtime();
    Trial rez2 = getSequential(trials, in_func, eps, r);
    double time4 = MPI_Wtime();
    std::cout << "time_p " << time2 - time1 << std::endl;
    std::cout << "time_s " << time4 - time3 << std::endl;
    EXPECT_NEAR(x_true, rez2.x, eps);
  }
}

TEST(Parallel_Operations_MPI, TestAGP_2) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<Trial> trials(2);
  int in_func = 1;
  double x_0 = -2, x_n = 3, x_true = 0;
  double eps = 0.0001, r = 2.5;
  if (rank == 0) {
    Trial f, s;
    f.x = x_0;
    f.z = Func(in_func, x_0);
    trials[0] = f;
    s.x = x_n;
    s.z = Func(in_func, x_n);
    trials[1] = s;
  }
  double time1 = MPI_Wtime();
  Trial rez1 = getParallelOperations(trials, in_func, eps, r);
  double time2 = MPI_Wtime();
  EXPECT_NEAR(x_true, rez1.x, eps);
  if (rank == 0) {
    double time3 = MPI_Wtime();
    Trial rez2 = getSequential(trials, in_func, eps, r);
    double time4 = MPI_Wtime();
    std::cout << "time_p " << time2 - time1 << std::endl;
    std::cout << "time_s " << time4 - time3 << std::endl;
    EXPECT_NEAR(x_true, rez2.x, eps);
  }
}

TEST(Parallel_Operations_MPI, TestAGP_3) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<Trial> trials(2);
  std::vector<int> coef = {2, 5, -1, 3, 0};
  int in_func = 2;
  double x_0 = 3, x_n = 7, x_true = 6.01286;
  double eps = 0.0001, r = 2.5;
  CoefIn(coef);
  if (rank == 0) {
    Trial f, s;
    f.x = x_0;
    f.z = Func(in_func, x_0);
    trials[0] = f;
    s.x = x_n;
    s.z = Func(in_func, x_n);
    trials[1] = s;
  }
  double time1 = MPI_Wtime();
  Trial rez1 = getParallelOperations(trials, in_func, eps, r);
  double time2 = MPI_Wtime();
  EXPECT_NEAR(x_true, rez1.x, eps);
  if (rank == 0) {
    double time3 = MPI_Wtime();
    Trial rez2 = getSequential(trials, in_func, eps, r);
    double time4 = MPI_Wtime();
    std::cout << "time_p " << time2 - time1 << std::endl;
    std::cout << "time_s " << time4 - time3 << std::endl;
    EXPECT_NEAR(x_true, rez2.x, eps);
  }
}

TEST(Parallel_Operations_MPI, TestAGP_4) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<Trial> trials(2);
  std::vector<int> coef = { 2, 1, 0, 0, 0 };
  int in_func = 2;
  double x_0 = -1, x_n = 5, x_true = 3*M_PI/2;
  double eps = 0.0001, r = 2.5;
  CoefIn(coef);
  if (rank == 0) {
    Trial f, s;
    f.x = x_0;
    f.z = Func(in_func, x_0);
    trials[0] = f;
    s.x = x_n;
    s.z = Func(in_func, x_n);
    trials[1] = s;
  }
  double time1 = MPI_Wtime();
  Trial rez1 = getParallelOperations(trials, in_func, eps, r);
  double time2 = MPI_Wtime();
  EXPECT_NEAR(x_true, rez1.x, eps);
  if (rank == 0) {
    double time3 = MPI_Wtime();
    Trial rez2 = getSequential(trials, in_func, eps, r);
    double time4 = MPI_Wtime();
    std::cout << "time_p " << time2 - time1 << std::endl;
    std::cout << "time_s " << time4 - time3 << std::endl;
    EXPECT_NEAR(x_true, rez2.x, eps);
  }
}

TEST(Parallel_Operations_MPI, TestAGP_5) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<Trial> trials(2);
  std::vector<int> coef = { -4, 1, 2, 3, 2 };
  int in_func = 2;
  double x_0 = 2, x_n = 4, x_true = 2.9064;
  double eps = 0.0001, r = 2.5;
  CoefIn(coef);
  if (rank == 0) {
    Trial f, s;
    f.x = x_0;
    f.z = Func(in_func, x_0);
    trials[0] = f;
    s.x = x_n;
    s.z = Func(in_func, x_n);
    trials[1] = s;
  }
  double time1 = MPI_Wtime();
  Trial rez1 = getParallelOperations(trials, in_func, eps, r);
  double time2 = MPI_Wtime();
  EXPECT_NEAR(x_true, rez1.x, eps);
  if (rank == 0) {
    double time3 = MPI_Wtime();
    Trial rez2 = getSequential(trials, in_func, eps, r);
    double time4 = MPI_Wtime();
    std::cout << "time_p " << time2 - time1 << std::endl;
    std::cout << "time_s " << time4 - time3 << std::endl;
    EXPECT_NEAR(x_true, rez2.x, eps);
  }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
