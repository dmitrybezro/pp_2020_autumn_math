// Copyright 2020 Vodeneev Mikhail
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./ops_mpi.h"

TEST(Parallel_Operations_MPI, Test_1) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n;
    if (rank == 0)
        n = 5;
    Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    ASSERT_EQ(n, 5);
}

TEST(Parallel_Operations_MPI, Test_2) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec(5);
    if (rank == 0) {
        for (int i = 0; i < 5; i++) {
            vec[i] = i;
        }
    }
    Bcast(vec.data(), 5, MPI_INT, 0, MPI_COMM_WORLD);
        ASSERT_EQ(vec[3], 3);
}


TEST(Parallel_Operations_MPI, Test_3) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> vec(1000);
    if (rank == 0) {
        for (int i = 0; i < 1000; i++) {
            vec[i] = i;
        }
    }
    Bcast(vec.data(), 1000, MPI_INT, 0, MPI_COMM_WORLD);
    ASSERT_EQ(vec[20], 20);
}

TEST(Parallel_Operations_MPI, Test_4) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> vec(1000000);
    if (rank == 0) {
        for (int i = 0; i < 1000000; i++) {
            vec[i] = i;
        }
    }
    Bcast(vec.data(), 1000000, MPI_INT, 0, MPI_COMM_WORLD);
    ASSERT_EQ(vec[999], 999);
}

TEST(Parallel_Operations_MPI, Test_5) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<float> vec(1000000);
    if (rank == 0) {
        for (int i = 0; i < 1000000; i++) {
            vec[i] = i;
        }
    }
    Bcast(vec.data(), 1000000, MPI_INT, 0, MPI_COMM_WORLD);
    ASSERT_EQ(vec[999], 999);
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
