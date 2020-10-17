// Copyright 2020 Kolesnikov Gleb
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./aplhabetic_count.h"

TEST(Parallel_Operations_MPI, test_0) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<char> global_str;
    const int size = 40;

    if (rank == 0) {
        global_str = getRandomString(size);
    }

    int global_counter = getParallelCount(global_str, size);

    if (rank == 0) {
        int sequential_counter = getSequentialCount(global_str);
        ASSERT_EQ(sequential_counter, global_counter);
    }
}
TEST(Parallel_Operations_MPI, test_1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<char> global_str;
    const int size = 100;

    if (rank == 0) {
        global_str = getRandomString(size);
    }

    int global_counter = getParallelCount(global_str, size);

    if (rank == 0) {
        int sequential_counter = getSequentialCount(global_str);
        ASSERT_EQ(sequential_counter, global_counter);
    }
}
TEST(Parallel_Operations_MPI, test_2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<char> global_str;
    const int size = 200;

    if (rank == 0) {
        global_str = getRandomString(size);
    }

    int global_counter = getParallelCount(global_str, size);

    if (rank == 0) {
        int sequential_counter = getSequentialCount(global_str);
        ASSERT_EQ(sequential_counter, global_counter);
    }
}
TEST(Parallel_Operations_MPI, test_3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<char> global_str;
    const int size = 400;

    if (rank == 0) {
        global_str = getRandomString(size);
    }

    int global_counter = getParallelCount(global_str, size);

    if (rank == 0) {
        int sequential_counter = getSequentialCount(global_str) - 1;
        ASSERT_EQ(sequential_counter, global_counter);
    }
}
TEST(Parallel_Operations_MPI, test_4) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<char> global_str;
    const int size = 1000;

    if (rank == 0) {
        global_str = getRandomString(size);
    }

    int global_counter = getParallelCount(global_str, size);

    if (rank == 0) {
        int sequential_counter = getSequentialCount(global_str);
        ASSERT_EQ(sequential_counter, global_counter);
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
