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
    matrix = getRandomMat(count_row, count_str);
    b = getRandomRes(count_str);
}
  try {
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
  }
  catch (char* war) {
    std::cerr << "war1: " << war << std::endl;
    ASSERT_ANY_THROW(getParallelMethod(matrix, b, count_row, count_str));
  }
  if (rank == 0) {
    try {
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
  catch (char* war) {
    std::cerr << "war2: " << war << std::endl;
    ASSERT_ANY_THROW(getSequentialMethod(matrix, b, count_row, count_str));
  }
  }
  }

TEST(Parallel_Operations_MPI, Test_MethodGaussa2) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int count_row = 6, count_str = 6;
  std::vector<double> matrix(count_row*count_str, 0), b(count_str, 0), x1(count_row, 0), x2(count_row, 0);
  double time1, time2, time3, time4, t_res1, t_res2, error1, error2;
  if (rank == 0) {
    matrix = getRandomMat(count_row, count_str);
    b = getRandomRes(count_str);
  }
  try {
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
  }
  catch (char* war) {
    std::cerr << "war1: " << war << std::endl;
    ASSERT_ANY_THROW(getParallelMethod(matrix, b, count_row, count_str));
  }
  if (rank == 0) {
    try {
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
    catch (char* war) {
      std::cerr << "war2: " << war << std::endl;
      ASSERT_ANY_THROW(getSequentialMethod(matrix, b, count_row, count_str));
    }
  }
  }

  TEST(Parallel_Operations_MPI, Test_MethodGaussa3) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int count_row = 13, count_str = 13;
  std::vector<double> matrix(count_row*count_str, 0), b(count_str, 0), x1(count_row, 0), x2(count_row, 0);
  double time1, time2, time3, time4, t_res1, t_res2, error1, error2;
  if (rank == 0) {
    matrix = getRandomMat(count_row, count_str);
  b = getRandomRes(count_str);
  }
  try {
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
  }
  catch (char* war) {
    std::cerr << "war1: " << war << std::endl;
    ASSERT_ANY_THROW(getParallelMethod(matrix, b, count_row, count_str));
  }
  if (rank == 0) {
    try {
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
    catch (char* war) {
      std::cerr << "war2: " << war << std::endl;
      ASSERT_ANY_THROW(getSequentialMethod(matrix, b, count_row, count_str));
    }
  }
  }

  TEST(Parallel_Operations_MPI, Test_MethodGaussa4) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int count_row = 20, count_str = 20;
  std::vector<double> matrix(count_row*count_str, 0), b(count_str, 0), x1(count_row, 0), x2(count_row, 0);
  double time1, time2, time3, time4, t_res1, t_res2, error1, error2;
  if (rank == 0) {
    matrix = getRandomMat(count_row, count_str);
    b = getRandomRes(count_str);
  }
  try {
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
  }
  catch (char* war) {
    std::cerr << "war1: " << war << std::endl;
    ASSERT_ANY_THROW(getParallelMethod(matrix, b, count_row, count_str));
  }
  if (rank == 0) {
    try {
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
    catch (char* war) {
      std::cerr << "war2: " << war << std::endl;
      ASSERT_ANY_THROW(getSequentialMethod(matrix, b, count_row, count_str));
    }
  }
  }

  TEST(Parallel_Operations_MPI, Test_MethodGaussa5) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int count_row = 25, count_str = 25;
  std::vector<double> matrix(count_row*count_str, 0), b(count_str, 0), x1(count_row, 0), x2(count_row, 0);
  double time1, time2, time3, time4, t_res1, t_res2, error1, error2;
  if (rank == 0) {
    matrix = getRandomMat(count_row, count_str);
    b = getRandomRes(count_str);
  }
  try {
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
  }
  catch (char* war) {
    std::cerr << "war1: " << war << std::endl;
    ASSERT_ANY_THROW(getParallelMethod(matrix, b, count_row, count_str));
  }
  if (rank == 0) {
    try {
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
    catch (char* war) {
      std::cerr << "war2: " << war << std::endl;
      ASSERT_ANY_THROW(getSequentialMethod(matrix, b, count_row, count_str));
    }
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
