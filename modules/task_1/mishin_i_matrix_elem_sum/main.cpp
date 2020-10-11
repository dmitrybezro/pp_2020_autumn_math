// Copyright 2020 Mishin Ilya
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./matrix_elem_sum.h"

TEST(Parallel_Operations_MPI, Test_Sum) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_matrix1D;
    const int count_size_cols = 10;
    const int count_size_rows = 10;

    if (rank == 0) {
        global_matrix1D = getRandomVector(count_size_rows * count_size_cols);
    }

    int global_sum = getParallelOperations(global_matrix1D, count_size_rows, count_size_cols , "+");

    if (rank == 0) {
        int reference_sum = getSequentialOperations(global_matrix1D, "+");
        ASSERT_EQ(reference_sum, global_sum);
    }
}

TEST(Parallel_Operations_MPI, Test_Diff) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_matrix1D;
    const int count_size_cols = 10;
    const int count_size_rows = 10;

    if (rank == 0) {
        global_matrix1D = getRandomVector(count_size_rows * count_size_cols);
    }

    int global_diff = getParallelOperations(global_matrix1D, count_size_rows, count_size_cols , "-");

    if (rank == 0) {
        int reference_diff = getSequentialOperations(global_matrix1D, "-");
        ASSERT_EQ(reference_diff, global_diff);
    }
}

TEST(Parallel_Operations_MPI, Test_Max) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_matrix1D;
    const int count_size_cols = 10;
    const int count_size_rows = 10;

    if (rank == 0) {
        global_matrix1D = getRandomVector(count_size_rows * count_size_cols);
    }

    int global_max = getParallelOperations(global_matrix1D, count_size_rows, count_size_cols , "max");

    if (rank == 0) {
        int reference_max = getSequentialOperations(global_matrix1D, "max");
        ASSERT_EQ(reference_max, global_max);
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
