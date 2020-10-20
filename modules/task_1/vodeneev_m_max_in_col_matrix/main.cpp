// Copyright 2020 Vodeneev Mikhail
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./ops_mpi.h"

TEST(Parallel_Operations_MPI, Test_Search_Max_50_50) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<std::vector<double>> test;
    int m = 50;  // row
    int n = 50;  // column
    if (rank == 0) {
        test = getRandomMatrix(m, n);
    }
    double time1 = MPI_Wtime();
    std::vector<double>res_p = getParallelOperations(test, m, n);
    double time2 = MPI_Wtime();
    double time_res1 = time2 - time1;
    if (rank == 0) {
        double time3 = MPI_Wtime();
        std::vector<double> res_s = getSeqOperations(test);
        double time4 = MPI_Wtime();
        double time_res2 = time4 - time3;
        std::cout << time_res1 << ' ' << time_res2;
        ASSERT_EQ(res_p, res_s);
    }
}

TEST(Parallel_Operations_MPI, Test_Search_Max_100_50) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<std::vector<double>> test;
    int m = 100;  // row
    int n = 50;  // column
    if (rank == 0) {
        test = getRandomMatrix(m, n);
    }
    double time1 = MPI_Wtime();
    std::vector<double>res_p = getParallelOperations(test, m, n);
    double time2 = MPI_Wtime();
    double time_res1 = time2 - time1;
    if (rank == 0) {
        double time3 = MPI_Wtime();
        std::vector<double> res_s = getSeqOperations(test);
        double time4 = MPI_Wtime();
        double time_res2 = time4 - time3;
        std::cout << time_res1 << ' ' << time_res2;
        ASSERT_EQ(res_p, res_s);
    }
}

TEST(Parallel_Operations_MPI, Test_Search_Max_999_999) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<std::vector<double>> test;
    int m = 999;  // row
    int n = 999;  // column
    if (rank == 0) {
        test = getRandomMatrix(m, n);
    }
    double time1 = MPI_Wtime();
    std::vector<double>res_p = getParallelOperations(test, m, n);
    double time2 = MPI_Wtime();
    double time_res1 = time2 - time1;
    if (rank == 0) {
        double time3 = MPI_Wtime();
        std::vector<double> res_s = getSeqOperations(test);
        double time4 = MPI_Wtime();
        double time_res2 = time4 - time3;
        std::cout << time_res1 << ' ' << time_res2;
        ASSERT_EQ(res_p, res_s);
    }
}

TEST(Parallel_Operations_MPI, Test_Search_Max_5000_5000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<std::vector<double>> test;
    int m = 5000;  // row
    int n = 5000;  // column
    if (rank == 0) {
        test = getRandomMatrix(m, n);
    }
    double time1 = MPI_Wtime();
    std::vector<double>res_p = getParallelOperations(test, m, n);
    double time2 = MPI_Wtime();
    double time_res1 = time2 - time1;
    if (rank == 0) {
        double time3 = MPI_Wtime();
        std::vector<double> res_s = getSeqOperations(test);
        double time4 = MPI_Wtime();
        double time_res2 = time4 - time3;
        std::cout << time_res1 << ' ' << time_res2;
        ASSERT_EQ(res_p, res_s);
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
