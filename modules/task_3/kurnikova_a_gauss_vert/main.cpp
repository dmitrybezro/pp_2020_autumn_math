// Copyright 2020 Kurnikova Anastasia
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "./ops_mpi.h"

TEST(Parallel_Operations_MPI, Test1) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    int n = 10;
    std::vector<double> a(n * n), res(n * n);
    for (int i = 0; i < n * n; i++)
        a[i] = 1;
    if (procrank == 0) {
        res = sequential(a, n, n);
        ASSERT_NO_THROW();
    }
}

TEST(Parallel_Operations_MPI, Test2) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    int n = 10;
    std::vector<double> a(n * n), res(n*n);
    if (procrank == 0) {
        for (int i = 0; i < n * n; i++)
            a[i] = 1;
    }
    res = parallel(a, n);
    if (procrank == 0) {
        ASSERT_NO_THROW();
    }
}

TEST(Parallel_Operations_MPI, Test3) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    int n = 10;
    std::vector<double> a(n * n), sres(n * n), pres(n * n);
    if (procrank == 0) {
        for (int i = 0; i < n * n; i++)
            a[i] = 0;
    }
    pres = parallel(a, n);
    if (procrank == 0) {
        sres = sequential(a, n, n);
        ASSERT_NO_THROW();
    }
}

TEST(Parallel_Operations_MPI, Test4) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    int n = 10;
    std::vector<double> a(n * n), sres(n * n), pres(n * n);
    if (procrank == 0) {
        for (int i = 0; i < n * n; i++)
            a[i] = 1;
    }
    pres = parallel(a, n);
    if (procrank == 0) {
        sres = sequential(a, n, n);
        ASSERT_NO_THROW();
    }
}

TEST(Parallel_Operations_MPI, Test5) {
    int procrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    int n = 1500;
    std::vector<double> a(n * n), sres(n * n), pres(n * n);
    if (procrank == 0) {
        for (int i = 0; i < n * n; i++)
            a[i] = 1;
    }
    double time1 = MPI_Wtime();
    pres = parallel(a, n);
    double time2 = MPI_Wtime();
    if (procrank == 0) {
        double res1 = time2 - time1;
        double time3 = MPI_Wtime();
        sres = sequential(a, n, n);
        double time4 = MPI_Wtime();
        double res2 = time4 - time3;
        std::cout << "Parallel time: " << res1 <<
             "\nSequential time: " << res2 << std::endl;
        ASSERT_NO_THROW();
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
