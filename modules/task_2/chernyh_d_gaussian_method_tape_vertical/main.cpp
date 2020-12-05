// Copyright 2018 Nesterov Alexander
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./method_gauss.h"

TEST(Parallel_Operations_MPI, Test_MethodGaussa1) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int count_row = 3, count_str = 3;
  std::vector<double> matrix(count_row*count_str, 0), b(count_str, 0), x1(count_row, 0), x2(count_row, 0);
  double time1, time2, time3, time4, t_res1, t_res2, error1, error2;
  if (rank == 0) {
    matrix = {
      3, 2, -5,
      2, -1, 3,
      1, 2, -1
    };
    b = {-1, 13, 9};
  }
  time1 = MPI_Wtime();
  x2 = getParallelMethod(matrix, b, count_row, count_str);
  time2 = MPI_Wtime();
  t_res1 = time2 - time1;
  error1 = SolutionCheck(matrix, b, x2, count_row, count_str);
  if (rank == 0) std::cout << "error1 " << error1 << std::endl;
  if (error1 <= FLT_EPSILON) error1 = 1;
  else
  error1 = 0;
  EXPECT_EQ(1, error1);
  if (rank == 0) {
    time3 = MPI_Wtime();
    x1 = getSequentialMethod(matrix, b, count_row, count_str);
    time4 = MPI_Wtime();
    t_res2 = time4 - time3;
    error2 = SolutionCheck(matrix, b, x1, count_row, count_str);
    std::cout << "error2 " << error2 << std::endl;
    std::cout << t_res1 << " " << t_res2 << std::endl;
    if (error2 <= FLT_EPSILON) error2 = 1;
    else
    error2 = 0;
    EXPECT_EQ(1, error2);
  }
}

TEST(Parallel_Operations_MPI, Test_MethodGaussa2) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int count_row = 3, count_str = 3;
  std::vector<double> matrix(count_row*count_str, 0), b(count_str, 0), x1(count_row, 0), x2(count_row, 0);
  double time1, time2, time3, time4, t_res1, t_res2, error1, error2;
  if (rank == 0) {
    matrix = {
      2, -3, 1,
      4, 3, -3,
      1, 7, -4
    };
    b = { 1, 3, 3 };
  }
  time1 = MPI_Wtime();
  x2 = getParallelMethod(matrix, b, count_row, count_str);
  time2 = MPI_Wtime();
  t_res1 = time2 - time1;
  error1 = SolutionCheck(matrix, b, x2, count_row, count_str);
  if (rank == 0) std::cout << "error1 " << error1 << std::endl;
  if (error1 <= FLT_EPSILON) error1 = 1;
  else
  error1 = 0;
  EXPECT_EQ(1, error1);
  if (rank == 0) {
    time3 = MPI_Wtime();
    x1 = getSequentialMethod(matrix, b, count_row, count_str);
    time4 = MPI_Wtime();
    t_res2 = time4 - time3;
    error2 = SolutionCheck(matrix, b, x1, count_row, count_str);
    std::cout << "error2 " << error2 << std::endl;
    std::cout << t_res1 << " " << t_res2 << std::endl;
    if (error2 <= FLT_EPSILON) error2 = 1;
    else
    error2 = 0;
    EXPECT_EQ(1, error2);
  }
}

TEST(Parallel_Operations_MPI, Test_MethodGaussa_Inconsistent_system) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int count_row = 3, count_str = 3;
  std::vector<double> matrix(count_row*count_str, 0), b(count_str, 0);
  if (rank == 0) {
    matrix = {
      1, 2, -4,
      2, 1, -5,
      1, -1, -1
    };
    b = { 1, -1, 3 };
  }
  if (rank == 0) {
    ASSERT_ANY_THROW(getSequentialMethod(matrix, b, count_row, count_str));
  }
}

TEST(Parallel_Operations_MPI, Test_MethodGaussa_Incorrect_dimensions) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int count_row = 3, count_str = 4;
  std::vector<double> matrix(count_row*count_str, 0), b(count_str, 0);
  ASSERT_ANY_THROW(getParallelMethod(matrix, b, count_row, count_str));
  if (rank == 0) {
    ASSERT_ANY_THROW(getSequentialMethod(matrix, b, count_row, count_str));
  }
}

TEST(Parallel_Operations_MPI, Test_MethodGaussa3) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int count_row = 4, count_str = 4;
  std::vector<double> matrix(count_row*count_str, 0), b(count_str, 0), x1(count_row, 0), x2(count_row, 0);
  double time1, time2, time3, time4, t_res1, t_res2, error1, error2;
  if (rank == 0) {
    matrix = {
      1, 1, 2, 3,
      1, 2, 3, -1,
      3, -1, -1, -2,
      2, 3, -1, -1
    };
    b = { 1, -4, -4, -6 };
  }
  time1 = MPI_Wtime();
  x2 = getParallelMethod(matrix, b, count_row, count_str);
  time2 = MPI_Wtime();
  t_res1 = time2 - time1;
  error1 = SolutionCheck(matrix, b, x2, count_row, count_str);
  if (rank == 0) std::cout << "error1 " << error1 << std::endl;
  if (error1 <= FLT_EPSILON) error1 = 1;
  else
  error1 = 0;
  EXPECT_EQ(1, error1);
  if (rank == 0) {
    time3 = MPI_Wtime();
    x1 = getSequentialMethod(matrix, b, count_row, count_str);
    time4 = MPI_Wtime();
    t_res2 = time4 - time3;
    error2 = SolutionCheck(matrix, b, x1, count_row, count_str);
    std::cout << "error2 " << error2 << std::endl;
    std::cout << t_res1 << " " << t_res2 << std::endl;
    if (error2 <= FLT_EPSILON) error2 = 1;
    else
    error2 = 0;
    EXPECT_EQ(1, error2);
  }
}

TEST(Parallel_Operations_MPI, Test_MethodGaussa4) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int count_row = 5, count_str = 5;
  std::vector<double> matrix(count_row*count_str, 0), b(count_str, 0), x1(count_row, 0), x2(count_row, 0);
  double time1, time2, time3, time4, t_res1, t_res2, error1, error2;
  if (rank == 0) {
    matrix = {
      3, -2, 1, 5, 6,
      5, 1, 3, 4, 5,
      1, 2, 3, 9, 8,
      7, 5, -5, 9, 1,
      1, -8, -2, 9, 1,
    };
    b = { 0, 0, 0, 0, 0};
  }
  time1 = MPI_Wtime();
  x2 = getParallelMethod(matrix, b, count_row, count_str);
  time2 = MPI_Wtime();
  t_res1 = time2 - time1;
  error1 = SolutionCheck(matrix, b, x2, count_row, count_str);
  if (rank == 0) std::cout << "error1 " << error1 << std::endl;
  if (error1 <= FLT_EPSILON) error1 = 1;
  else
  error1 = 0;
  EXPECT_EQ(1, error1);
  if (rank == 0) {
    time3 = MPI_Wtime();
    x1 = getSequentialMethod(matrix, b, count_row, count_str);
    time4 = MPI_Wtime();
    t_res2 = time4 - time3;
    error2 = SolutionCheck(matrix, b, x1, count_row, count_str);
    std::cout << "error2 " << error2 << std::endl;
    std::cout << t_res1 << " " << t_res2 << std::endl;
    if (error2 <= FLT_EPSILON) error2 = 1;
    else
    error2 = 0;
    EXPECT_EQ(1, error2);
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
