// Copyright 2020 Kudriavtsev Alexander
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./ops_mpi.h"

TEST(Parallel_Operations_MPI, Test_Seq_20) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int n = 20;
    const int NMax = 2000;
    const double eps = 0.00001;
    std::vector<double> A;
    std::vector<double> b;
    if (rank == 0) {
        A = getRandomMatrix(n);
        b = getRandomVector(n);
        std::vector<double> x = sequentialIterMethod(A, b, n, eps, NMax);
        double disp = discrepancyNorm(x, A, b);
        if (disp > 1.0) {
            std::cout << "Oh, ****, I'm sorry" << std::endl;
        }
        else ASSERT_EQ(true, disp < 1.0);
    }
}

TEST(Parallel_Operations_MPI, Test_Par_20) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int n = 20;
    const int NMax = 2000;
    const double eps = 0.00001;
    std::vector<double> A;
    std::vector<double> b;

    if (rank == 0) {
        A = getRandomMatrix(n);
        b = getRandomVector(n);
    }
    std::vector<double> x1 = parallelIterMethod(A, b, n, eps, NMax);
    if (rank == 0) {
        double disp = discrepancyNorm(x1, A, b);
        if (disp > 1.0) {
            std::cout << "Oh, ****, I'm sorry" << std::endl;
        }
        else ASSERT_EQ(true, disp < 1.0);
    }
}

TEST(Parallel_Operations_MPI, Test_Par_40) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int n = 40;
    const int NMax = 2000;
    const double eps = 0.00001;
    std::vector<double> A;
    std::vector<double> b;

    if (rank == 0) {
        A = getRandomMatrix(n);
        b = getRandomVector(n);
    }
    std::vector<double> x = parallelIterMethod(A, b, n, eps, NMax);
    if (rank == 0) {
        double disp = discrepancyNorm(x, A, b);
        if (disp > 1.0) {
            std::cout << "Oh, ****, I'm sorry" << std::endl;
        }
        else ASSERT_EQ(true, disp < 1.0);
    }
}

TEST(Parallel_Operations_MPI, Test_Par_100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int n = 100;
    const int NMax = 2000;
    const double eps = 0.00001;
    std::vector<double> A;
    std::vector<double> b;

    if (rank == 0) {
        A = getRandomMatrix(n);
        b = getRandomVector(n);
    }
    std::vector<double> x = parallelIterMethod(A, b, n, eps, NMax);
    if (rank == 0) {
        double disp = discrepancyNorm(x, A, b);
        if (disp > 1.0) {
            std::cout << "Oh, ****, I'm sorry" << std::endl;
        }
        else ASSERT_EQ(true, disp < 1.0);
    }
}

TEST(Parallel_Operations_MPI, Test_Seq_40) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int n = 40;
    const int NMax = 2000;
    const double eps = 0.00001;
    std::vector<double> A;
    std::vector<double> b;

    if (rank == 0) {
        A = getRandomMatrix(n);
        b = getRandomVector(n);
        std::vector<double> x = sequentialIterMethod(A, b, n, eps, NMax);
        double disp = discrepancyNorm(x, A, b);
        if (disp > 1.0) {
            std::cout << "Oh, ****, I'm sorry" << std::endl;
        }
        else ASSERT_EQ(true, disp < 1.0);
    }
}

TEST(Parallel_Operations_MPI, Test_Seq_100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int n = 100;
    const int NMax = 2000;
    const double eps = 0.00001;
    std::vector<double> A;
    std::vector<double> b;

    if (rank == 0) {
        A = getRandomMatrix(n);
        b = getRandomVector(n);
        std::vector<double> x = sequentialIterMethod(A, b, n, eps, NMax);
        double disp = discrepancyNorm(x, A, b);
        if (disp > 1.0) {
            std::cout << "Oh, ****, I'm sorry" << std::endl;
        }
        else ASSERT_EQ(true, disp < 1.0);
    }
}

TEST(Parallel_Operations_MPI, Test_Seq_and_Par_100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int n = 100;
    const int NMax = 2000;
    const double eps = 0.00001;
    std::vector<double> A;
    std::vector<double> b;

    if (rank == 0) {
        A = getRandomMatrix(n);
        b = getRandomVector(n);
    }
    std::vector<double> x1 = parallelIterMethod(A, b, n, eps, NMax);
    if (rank == 0) {
        std::vector<double> x2 = sequentialIterMethod(A, b, n, eps, NMax);
        double disp1 = discrepancyNorm(x1, A, b);
        double disp2 = discrepancyNorm(x1, A, b);
        if (disp1 >= 1.0 && disp2 >= 1.0) {
            std::cout << "Oh, ****, I'm sorry" << std::endl;
        }
        else {
            ASSERT_EQ(true, disp1 < 1);
            ASSERT_EQ(true, disp2 < 1);
        }
    }
}

TEST(Parallel_Operations_MPI, Test_Seq_and_Par_300) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int n = 300;
    const int NMax = 2000;
    const double eps = 0.00001;
    std::vector<double> A;
    std::vector<double> b;

    if (rank == 0) {
        A = getRandomMatrix(n);
        b = getRandomVector(n);
    }
    std::vector<double> x1 = parallelIterMethod(A, b, n, eps, NMax);
    if (rank == 0) {
        std::vector<double> x2 = sequentialIterMethod(A, b, n, eps, NMax);
        double disp1 = discrepancyNorm(x1, A, b);
        double disp2 = discrepancyNorm(x1, A, b);
        if (disp1 >= 1.0 && disp2 >= 1.0) {
            std::cout << "Oh, ****, I'm sorry" << std::endl;
        }
        else {
            ASSERT_EQ(true, disp1 < 1);
            ASSERT_EQ(true, disp2 < 1);
        }
    }
}

TEST(Parallel_Operations_MPI, Test_Seq_and_Par_727) {  // Simple Number
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int n = 727;
    const int NMax = 2000;
    const double eps = 0.00001;
    std::vector<double> A;
    std::vector<double> b;

    if (rank == 0) {
        A = getRandomMatrix(n);
        b = getRandomVector(n);
    }
    std::vector<double> x1 = parallelIterMethod(A, b, n, eps, NMax);
    if (rank == 0) {
        std::vector<double> x2 = sequentialIterMethod(A, b, n, eps, NMax);
        double disp1 = discrepancyNorm(x1, A, b);
        double disp2 = discrepancyNorm(x1, A, b);
        if (disp1 >= 1.0 && disp2 >= 1.0) {
            std::cout << "Oh, ****, I'm sorry" << std::endl;
        }
        else {
            ASSERT_EQ(true, disp1 < 1);
            ASSERT_EQ(true, disp2 < 1);
        }
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
