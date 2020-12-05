// Copyright 2020 Kurnikova Anastasia
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "./ops_mpi.h"

TEST(Parallel_Operations_MPI, Test1) {
    int procrank, done;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    int n = 3;
    std::vector<double> a(n * n);
    a[0] = 10;
    a[1] = 1;
    a[2] = -1;
    a[3] = 1;
    a[4] = 10;
    a[5] = -1;
    a[6] = -1;
    a[7] = 1;
    a[8] = 10;
    if (procrank == 0) {
        done = correctness(a, n);
        ASSERT_EQ(done, 1);
    }
}

TEST(Parallel_Operations_MPI, Test2) {
    int procrank, done;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    int n = 3;
    std::vector<double> a(n * n);
    a[0] = 5;
    a[1] = 4;
    a[2] = 3;
    a[3] = -1;
    a[4] = 0;
    a[5] = 2;
    a[6] = 1;
    a[7] = 1;
    a[8] = 0;
    if (procrank == 0) {
        done = correctness(a, n);
        ASSERT_EQ(done, 0);
    }
}

TEST(Parallel_Operations_MPI, Test3) {
    int procrank, done;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    int n = 3;
    const double e = 0.001;
    std::vector<double> a(n * n);
    a[0] = 10;
    a[1] = 1;
    a[2] = -1;
    a[3] = 1;
    a[4] = 10;
    a[5] = -1;
    a[6] = -1;
    a[7] = 1;
    a[8] = 10;
    std::vector<double> b(n);
    b[0] = 11;
    b[1] = 10;
    b[2] = 10;
    if (procrank == 0) {
        done = sequential(a, b, n, e);
        ASSERT_EQ(done, 1);
    }
}

TEST(Parallel_Operations_MPI, Test4) {
    int procrank, done;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    int n = 3;
    const double e = 0.001;
    std::vector<double> a(n * n);
    a[0] = 10;
    a[1] = 1;
    a[2] = -1;
    a[3] = 1;
    a[4] = 10;
    a[5] = -1;
    a[6] = -1;
    a[7] = 1;
    a[8] = 10;
    std::vector<double> b(n);
    b[0] = 11;
    b[1] = 10;
    b[2] = 10;
    done = parallel(a, b, n, e);
    if (procrank == 0) {
        ASSERT_EQ(done, 1);
    }
}

TEST(Parallel_Operations_MPI, Test5) {
    int procrank, done;
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    int n = 3;
    const double e = 0.001;
    std::vector<double> a(n * n);
    a[0] = 10;
    a[1] = 1;
    a[2] = -1;
    a[3] = 1;
    a[4] = 10;
    a[5] = -1;
    a[6] = -1;
    a[7] = 1;
    a[8] = 10;
    std::vector<double> b(n);
    b[0] = 11;
    b[1] = 10;
    b[2] = 10;
    double time1 = MPI_Wtime();
    done = parallel(a, b, n, e);
    double time2 = MPI_Wtime();
    double res2 = time2 - time1;
    if (procrank == 0) {
        double time3 = MPI_Wtime();
        done = sequential(a, b, n, e);
        double time4 = MPI_Wtime();
        double res1 = time4 - time3;
        std::cout << res1 << ' ' << res2 << std::endl;
        ASSERT_EQ(done, 1);
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
