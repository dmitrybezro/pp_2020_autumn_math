// Copyright 2020 Kurnikova Anastasia
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include "./ops_mpi.h"

TEST(Parallel_Operations_MPI, Test1) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    double parans = getParallelOperations(f, 0, 5, 13);
    if (procrank == 0) {
        double seqans = getSequentialOperations(f, 0, 5, 13);
        ASSERT_EQ(seqans, parans);
    }
}

TEST(Parallel_Operations_MPI, Test2) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    double parans = getParallelOperations(f, -4, 0, 2);
    if (procrank == 0) {
        double seqans = getSequentialOperations(f, -4, 0, 2);
        ASSERT_EQ(seqans, parans);
    }
}

TEST(Parallel_Operations_MPI, Test3) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    double parans = getParallelOperations(f, 12, 12, 12);
    if (procrank == 0) {
        double seqans = getSequentialOperations(f, 12, 12, 12);
        ASSERT_EQ(seqans, parans);
    }
}

TEST(Parallel_Operations_MPI, Test4) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    double parans = getParallelOperations(f, 5, 20, 10);
    if (procrank == 0) {
        double seqans = getSequentialOperations(f, 5, 20, 10);
        ASSERT_EQ(seqans, parans);
    }
}

TEST(Parallel_Operations_MPI, Test5) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    double parans = getParallelOperations(f, 11, 15, 3);
    if (procrank == 0) {
        double seqans = getSequentialOperations(f, 11, 15, 3);
        ASSERT_EQ(seqans, parans);
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
