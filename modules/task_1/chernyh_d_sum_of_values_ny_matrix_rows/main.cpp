// Copyright 2018 Nesterov Alexander
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./matrix_sum.h"

TEST(Parallel_Operations_MPI, Test_Sum) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<int> global_mat;
  const int count_row = 5, count_str = 5;
  if (rank == 0) {
    global_mat = getRandomMat(count_row, count_str);
  }
  double time1 = MPI_Wtime();
  std::vector<int> global_sum = getParallelOperations(global_mat, count_row, count_str);
  double time2 = MPI_Wtime();
  double t_res1 = time2 - time1;
  if (rank == 0) {
    double time3 = MPI_Wtime();
    std::vector<int> reference_sum = getSequentialOperations(global_mat, count_row, count_str);
    double time4 = MPI_Wtime();
    double t_res2 = time4 - time3;
    std::cout << t_res1 << " " << t_res2 << std::endl;
    ASSERT_EQ(reference_sum, global_sum);
  }
}

TEST(Parallel_Operations_MPI, Test_Sum_2) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<int> global_mat;
  const int count_row = 100, count_str = 100;
  if (rank == 0) {
    global_mat = getRandomMat(count_row, count_str);
  }
  double time1 = MPI_Wtime();
  std::vector<int> global_sum = getParallelOperations(global_mat, count_row, count_str);
  double time2 = MPI_Wtime();
  double t_res1 = time2 - time1;
  if (rank == 0) {
    double time3 = MPI_Wtime();
    std::vector<int> reference_sum = getSequentialOperations(global_mat, count_row, count_str);
    double time4 = MPI_Wtime();
    double t_res2 = time4 - time3;
    std::cout << t_res1 << " " << t_res2 << std::endl;
    ASSERT_EQ(reference_sum, global_sum);
}
}

TEST(Parallel_Operations_MPI, Test_Sum_3) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<int> global_mat;
  const int count_row = 105, count_str = 213;
  if (rank == 0) {
    global_mat = getRandomMat(count_row, count_str); }
  double time1 = MPI_Wtime();
  std::vector<int> global_sum = getParallelOperations(global_mat, count_row, count_str);
  double time2 = MPI_Wtime();
  double t_res1 = time2 - time1;
  if (rank == 0) {
    double time3 = MPI_Wtime();
    std::vector<int> reference_sum = getSequentialOperations(global_mat, count_row, count_str);
    double time4 = MPI_Wtime();
    double t_res2 = time4 - time3;
    std::cout << t_res1 << " " << t_res2 << std::endl;
    ASSERT_EQ(reference_sum, global_sum); }
}

TEST(Parallel_Operations_MPI, Test_Sum_4) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<int> global_mat;
  const int count_row = 1002, count_str = 2001;
  if (rank == 0) {
    global_mat = getRandomMat(count_row, count_str);
  }
  double time1 = MPI_Wtime();
  std::vector<int> global_sum = getParallelOperations(global_mat, count_row, count_str);
  double time2 = MPI_Wtime();
  double t_res1 = time2 - time1;
  if (rank == 0) {
    double time3 = MPI_Wtime();
    std::vector<int> reference_sum = getSequentialOperations(global_mat, count_row, count_str);
    double time4 = MPI_Wtime();
    double t_res2 = time4 - time3;
    std::cout << t_res1 << " " << t_res2 << std::endl;
    ASSERT_EQ(reference_sum, global_sum);
  }
}

TEST(Parallel_Operations_MPI, Test_Sum_5) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<int> global_mat;
  const int count_row = 10000, count_str = 10000;
  if (rank == 0) {
    global_mat = getRandomMat(count_row, count_str);
  }
  double time1 = MPI_Wtime();
  std::vector<int> global_sum = getParallelOperations(global_mat, count_row, count_str);
  double time2 = MPI_Wtime();
  double t_res1 = time2 - time1;
  if (rank == 0) {
    double time3 = MPI_Wtime();
    std::vector<int> reference_sum = getSequentialOperations(global_mat, count_row, count_str);
    double time4 = MPI_Wtime();
    double t_res2 = time4 - time3;
    std::cout << t_res1 << " " << t_res2 << std::endl;
    ASSERT_EQ(reference_sum, global_sum);
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
