// Copyright 2018 Smirnov Alexander
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <tuple>
#include "./sparse_matrix_multiplication_ccs.h"

TEST(Parallel_Operations_MPI, Test_Transpose) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> Avalues = {3.0, 7.0, 8.0, 9.0, 15.0, 16.0};
    std::vector<int> Acols = {1, 3, 2, 0, 2, 3};
    std::vector<int> Apointers = {0, 2, 3, 3, 6};

    std::vector<double> Bvalues = {9.0, 3.0, 8.0, 15.0, 7.0, 16.0};
    std::vector<int> Bcols = {3, 0, 1, 3, 0, 3};
    std::vector<int> Bpointers = {0, 1, 2, 4, 6};

    if (rank == 0) {
        auto newA = transpose(Avalues, Acols, Apointers);
        ASSERT_EQ(std::get<0>(newA), Bvalues);
        ASSERT_EQ(std::get<1>(newA), Bcols);
        ASSERT_EQ(std::get<2>(newA), Bpointers);
    }
}

TEST(Parallel_Operations_MPI, Test_DoubleTranspose) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> Avalues = { 3.0, 7.0, 8.0, 9.0, 15.0, 16.0 };
    std::vector<int> Acols = { 1, 3, 2, 0, 2, 3 };
    std::vector<int> Apointers = { 0, 2, 3, 3, 6 };

    if (rank == 0) {
        auto newA = transpose(Avalues, Acols, Apointers);
        auto newAA = transpose(std::get<0>(newA), std::get<1>(newA), std::get<2>(newA));
        ASSERT_EQ(std::get<0>(newAA), Avalues);
        ASSERT_EQ(std::get<1>(newAA), Acols);
        ASSERT_EQ(std::get<2>(newAA), Apointers);
    }
}

TEST(Parallel_Operations_MPI, Test_TransposeIdentityMatrix) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> Avalues = { 1.0, 1.0, 1.0, 1.0 };
    std::vector<int> Acols = { 0, 1, 2, 3 };
    std::vector<int> Apointers = { 0, 1, 2, 3, 4 };

    if (rank == 0) {
        auto newA = transpose(Avalues, Acols, Apointers);
        ASSERT_EQ(std::get<0>(newA), Avalues);
        ASSERT_EQ(std::get<1>(newA), Acols);
        ASSERT_EQ(std::get<2>(newA), Apointers);
    }
}

TEST(Parallel_Operations_MPI, Test_SequentialMultiplication) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> Avalues = { 3.0, 7.0, 8.0, 9.0, 15.0, 16.0 };
    std::vector<int> Acols = { 1, 3, 2, 0, 2, 3 };
    std::vector<int> Apointers = { 0, 2, 3, 3, 6 };

    std::vector<double> Bvalues = { 9.0, 3.0, 8.0, 15.0, 7.0, 16.0 };
    std::vector<int> Bcols = { 3, 0, 1, 3, 0, 3 };
    std::vector<int> Bpointers = { 0, 1, 2, 4, 6 };

    std::vector<double> Cvalues = { 63.0, 129.0, 112.0, 144.0, 27.0, 240.0, 319.0 };
    std::vector<int> Ccols = { 0, 2, 3, 0, 1, 2, 3 };
    std::vector<int> Cpointers = { 0, 3, 3, 3, 7 };

    if (rank == 0) {
        auto newC = getSequentialOperations(Avalues, Acols, Apointers, Bvalues, Bcols, Bpointers);
        ASSERT_EQ(std::get<0>(newC), Cvalues);
        ASSERT_EQ(std::get<1>(newC), Ccols);
        ASSERT_EQ(std::get<2>(newC), Cpointers);
    }
}

TEST(Parallel_Operations_MPI, Test_ParallelMultiplication) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 4;
    std::vector<double> Avalues = { 3.0, 7.0, 8.0, 9.0, 15.0, 16.0 };
    std::vector<int> Acols = { 1, 3, 2, 0, 2, 3 };
    std::vector<int> Apointers = { 0, 2, 3, 3, 6 };

    std::vector<double> Bvalues = { 9.0, 3.0, 8.0, 15.0, 7.0, 16.0 };
    std::vector<int> Bcols = { 3, 0, 1, 3, 0, 3 };
    std::vector<int> Bpointers = { 0, 1, 2, 4, 6 };

    std::vector<double> Cvalues = { 63.0, 129.0, 112.0, 144.0, 27.0, 240.0, 319.0 };
    std::vector<int> Ccols = { 0, 2, 3, 0, 1, 2, 3 };
    std::vector<int> Cpointers = { 0, 3, 3, 3, 7 };

    auto newC = getParallelOperations(size, Avalues, Acols, Apointers, Bvalues, Bcols, Bpointers);

    if (rank == 0) {
        ASSERT_EQ(std::get<0>(newC), Cvalues);
        ASSERT_EQ(std::get<1>(newC), Ccols);
        ASSERT_EQ(std::get<2>(newC), Cpointers);
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
