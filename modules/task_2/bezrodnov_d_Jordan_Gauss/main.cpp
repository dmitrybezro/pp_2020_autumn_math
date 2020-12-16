// Copyright 2020 Bezrodnov Dmitry
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include<mpi.h>
#include<cmath>
#include <iostream>
#include <vector>
#include "./Jordan_Gauss.h"

TEST(Parallel_Operations_MPI, Test_1_Matrix_3_3) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    int size = 3;
    std::vector<double> MainMatr = { 3, 2, 1, 2, -1, 7, -1, 5, -1, 4, 23, 5 };
    std::vector<double> X_parall = ParallelJordanGauss(MainMatr, size);
    if (RANK == 0) {
        std::vector<double> X_inv = MultiInverseMatrix(MainMatr, size);
        for (int i = 0; i < X_parall.size(); i++) {
           ASSERT_NEAR(X_parall[i], X_inv[i], 0.01);
        }
    }
}

TEST(Parallel_Operations_MPI, Test_2_Matrix_4_4) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    int size = 4;
    std::vector<double> MainMatr = { 3, 3, 2, 1, 3, 1, 1, 3, 6, 5, 4, 3, 3, 1, 2, 2, 6, 2, 1, 6 };
    std::vector<double> X_parall = ParallelJordanGauss(MainMatr, size);
    if (RANK == 0) {
        std::vector<double> X_inv = MultiInverseMatrix(MainMatr, size);
        for (int i = 0; i < X_parall.size(); i++) {
           ASSERT_NEAR(X_parall[i], X_inv[i], 0.01);
        }
    }
}

TEST(Parallel_Operations_MPI, Test_3_Matrix_5_5) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    int size = 5;
    std::vector<double> MainMatr = { 1, 2, 2, 2, 7,
                                2, -1, 4, 1, 2,
                                2, 1, -3, 2, 3,
                                -3, 2, 5, -1, -1,
                                -2, -1, 4, -2, 5,
                                8, 3, 6, 19, 0 };
    std::vector<double> X_parall = ParallelJordanGauss(MainMatr, size);
    if (RANK == 0) {
        std::vector<double> X_inv = MultiInverseMatrix(MainMatr, size);
        for (int i = 0; i < X_parall.size(); i++) {
           ASSERT_NEAR(X_parall[i], X_inv[i], 0.01);
        }
    }
}

TEST(Parallel_Operations_MPI, Test_4_Matrix_7_7) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    int size = 7;
    std::vector<double> MainMatr = getRandomMatrix(size);
    std::vector<double> X_parall = ParallelJordanGauss(MainMatr, size);
    if (RANK == 0) {
        std::vector<double> X_sequen = SequenJordanGauss(MainMatr, size);
        std::vector<double> X_inv = MultiInverseMatrix(MainMatr, size);
        for (int i = 0; i < X_parall.size(); i++) {
           ASSERT_NEAR(X_parall[i], X_inv[i], 0.01);
        }
    }
}

TEST(Parallel_Operations_MPI, Test_5_Matrix_8_8) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    int size = 8;
    std::vector<double> MainMatr = getRandomMatrix(size);
    std::vector<double> X_parall = ParallelJordanGauss(MainMatr, size);
    if (RANK == 0) {
        std::vector<double> X_inv = MultiInverseMatrix(MainMatr, size);
        for (int i = 0; i < X_parall.size(); i++) {
           ASSERT_NEAR(X_parall[i], X_inv[i], 0.01);
        }
    }
}

/*TEST(Parallel_Operations_MPI, Test_5_Matrix_100_100) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    int size = 2000;
    std::vector<double> MainMatr = getRandomMatrix(size);
    double start_1 = MPI_Wtime();
    std::vector<double> X_parall = ParallelJordanGauss(MainMatr, size);
    double finish_1 = MPI_Wtime();
    if (RANK == 0) {
        double start_2 = MPI_Wtime();
        std::vector<double> X_sequen = SequenJordanGauss(MainMatr, size);
        double finish_2 = MPI_Wtime();
        std::cout << "Time parall = " << finish_1 - start_1 <<std::endl;
        std::cout << "Time sequen = " << finish_2 - start_2 <<std::endl;
        for (int i = 0; i < X_parall.size(); i++) {
           ASSERT_NEAR(X_parall[i], X_sequen[i], 0.01);
        }
    }
}*/

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
