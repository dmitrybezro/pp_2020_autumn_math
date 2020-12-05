// Copyright 2020 Bezrodnov Dmitry
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./min_elem_matrix.h"

TEST(Parallel_Operations_MPI, Matrix_5_5) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    const int rows = 5;
    const int cols = 5;
    std::vector<int> matr;
    if (RANK == 0) {
        matr = getRandomMatrix(rows, cols);
    }

    int min_parall = getParallelOperations(matr, rows, cols);

    if (RANK == 0) {
        int min_sequent = getSequentialOperations(matr);
        ASSERT_EQ(min_sequent, min_parall);
    }
}

TEST(Parallel_Operations_MPI, Test_10_20) {
     int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    const int rows = 10;
    const int cols = 20;
    std::vector<int> matr;
    if (RANK == 0) {
        matr = getRandomMatrix(rows, cols);
    }

    int min_parall = getParallelOperations(matr, rows, cols);

    if (RANK == 0) {
        int min_sequent = getSequentialOperations(matr);
        ASSERT_EQ(min_sequent, min_parall);
    }
}

TEST(Parallel_Operations_MPI, Test_101_202) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    const int rows = 101;
    const int cols = 202;
    std::vector<int> matr;
    if (RANK == 0) {
        matr = getRandomMatrix(rows, cols);
    }

    int min_parall = getParallelOperations(matr, rows, cols);

    if (RANK == 0) {
        int min_sequent = getSequentialOperations(matr);
        ASSERT_EQ(min_sequent, min_parall);
    }
}

TEST(Parallel_Operations_MPI, Test_1000_2000) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    const int rows = 1000;
    const int cols = 2000;
    std::vector<int> matr;
    if (RANK == 0) {
        matr = getRandomMatrix(rows, cols);
    }

    int min_parall = getParallelOperations(matr, rows, cols);

    if (RANK == 0) {
        int min_sequent = getSequentialOperations(matr);
        ASSERT_EQ(min_sequent, min_parall);
    }
}

TEST(Parallel_Operations_MPI, Test_2002_3001) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    const int rows = 2002;
    const int cols = 3001;
    std::vector<int> matr;
    if (RANK == 0) {
        matr = getRandomMatrix(rows, cols);
    }

    int min_parall = getParallelOperations(matr, rows, cols);

    if (RANK == 0) {
        int min_sequent = getSequentialOperations(matr);
        ASSERT_EQ(min_sequent, min_parall);
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

