// Copyright 2020 Lebedev Andrew
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./matrix_max_elem.h"

TEST(Parallel_Operations_MPI, Test_0_0) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_mtx;
    const int rows = 0;
    const int cols = 0;
    double time_start, time_stop;
    if (rank == 0) {
        global_mtx = getRandomMatrix(rows * cols);
    }

    int parallel_max = getParallelOperations(global_mtx, rows * cols);

    if (rank == 0) {
        time_start = MPI_Wtime();
        int sequential_max = getSequentialOperations(global_mtx);
        time_stop = MPI_Wtime();
        printf("Sequential time: %3.20f\n", time_stop - time_start);
        ASSERT_EQ(sequential_max, parallel_max);
    }
}

TEST(Parallel_Operations_MPI, Test_10_10) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_mtx;
    const int rows = 10;
    const int cols = 10;
    double time_start, time_stop;
    if (rank == 0) {
        global_mtx = getRandomMatrix(rows * cols);
    }

    int parallel_max = getParallelOperations(global_mtx, rows * cols);

    if (rank == 0) {
        time_start = MPI_Wtime();
        int sequential_max = getSequentialOperations(global_mtx);
        time_stop = MPI_Wtime();
        printf("Sequential time: %3.20f\n", time_stop - time_start);
        ASSERT_EQ(sequential_max, parallel_max);
    }
}

TEST(Parallel_Operations_MPI, Test_20_10) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_mtx;
    const int rows = 10;
    const int cols = 10;
    double time_start, time_stop;
    if (rank == 0) {
        global_mtx = getRandomMatrix(rows * cols);
    }

    int parallel_max = getParallelOperations(global_mtx, rows * cols);

    if (rank == 0) {
        time_start = MPI_Wtime();
        int sequential_max = getSequentialOperations(global_mtx);
        time_stop = MPI_Wtime();
        printf("Sequential time: %3.20f\n", time_stop - time_start);
        ASSERT_EQ(sequential_max, parallel_max);
    }
}

TEST(Parallel_Operations_MPI, Test_100_100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_mtx;
    const int rows = 100;
    const int cols = 100;
    double time_start, time_stop;
    if (rank == 0) {
        global_mtx = getRandomMatrix(rows * cols);
    }

    int parallel_max = getParallelOperations(global_mtx, rows * cols);

    if (rank == 0) {
        time_start = MPI_Wtime();
        int sequential_max = getSequentialOperations(global_mtx);
        time_stop = MPI_Wtime();
        printf("Sequential time: %3.20f\n", time_stop - time_start);
        ASSERT_EQ(sequential_max, parallel_max);
    }
}

TEST(Parallel_Operations_MPI, Test_1000_1000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_mtx;
    const int rows = 1000;
    const int cols = 1000;
    double time_start, time_stop;
    if (rank == 0) {
        global_mtx = getRandomMatrix(rows * cols);
    }

    int parallel_max = getParallelOperations(global_mtx, rows * cols);

    if (rank == 0) {
        time_start = MPI_Wtime();
        int sequential_max = getSequentialOperations(global_mtx);
        time_stop = MPI_Wtime();
        printf("Sequential time: %3.20f\n", time_stop - time_start);
        ASSERT_EQ(sequential_max, parallel_max);
    }
}

TEST(Parallel_Operations_MPI, Test_1000_10000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_mtx;
    const int rows = 1000;
    const int cols = 10000;
    double time_start, time_stop;
    if (rank == 0) {
        global_mtx = getRandomMatrix(rows * cols);
    }

    int parallel_max = getParallelOperations(global_mtx, rows * cols);

    if (rank == 0) {
        time_start = MPI_Wtime();
        int sequential_max = getSequentialOperations(global_mtx);
        time_stop = MPI_Wtime();
        printf("Sequential time: %3.20f\n", time_stop - time_start);
        ASSERT_EQ(sequential_max, parallel_max);
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
