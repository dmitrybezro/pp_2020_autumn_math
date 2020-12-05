// Copyright 2020 Kolesnikov Gleb
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <ctime>
#include <iostream>
#include <cmath>
#include "./ops_mpi.h"
double f4(std::vector<double> coords) {
    return (coords[0] + coords[1] + coords[2]);
}
double f3(std::vector<double> coords) {
    return (coords[0] * sin(coords[1]));
}
double f2(std::vector<double> coords) {
    return (coords[0] * coords[1] * coords[2]);
}
double f1(std::vector<double> coords) {
    return (1);
}
double f0(std::vector<double> coords) {
    return (0);
}
TEST(Trapezoid_Rule, Test_f) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double error = 0.01;
    std::vector<double> x(3);
    std::vector<double> y(3);
    x = { 0, 0, 0 };
    y = { 2, 2, 2 };
    const int n = 50;
    double res = trapezoidParallelRule(f4, x, y, n);
    if (rank == 0) {
        ASSERT_EQ(abs(res - 24.0) < error, 1);
    }
}
TEST(Trapezoid_Rule, Test_f3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int dimension = 3;
    double error = 0.001;
    int t;
    std::vector<double> x(dimension);
    std::vector<double> y(dimension);
    x = { 0, 0, 0 };
    y = { 3, 3, 3 };
    const int n = 150;
    t = clock();
    double res = trapezoidParallelRule(f3, x, y, n);
    if (rank == 0) {
        std::cout << "parallel res: " << res << std::endl;
        std::cout << "parallel time: " << clock() - t << std::endl;
        ASSERT_EQ(abs(res - 26.86489) < error, 1);
    }
}
TEST(Trapezoid_Rule, Test_f2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int dimension = 3;
    double error = 0.001;
    int t;
    std::vector<double> x(dimension);
    std::vector<double> y(dimension);
    x = { 0, 0, 0 };
    y = { 2, 2, 2 };
    const int n = 50;
    t = clock();
    double res = trapezoidParallelRule(f2, x, y, n);
    if (rank == 0) {
        std::cout << "parallel res: " << res << std::endl;
        std::cout << "parallel time: " << clock() - t << std::endl;
        t = clock();
        ASSERT_EQ(abs(res - 8.0) < error, 1);
        std::cout << "seq res: " << SequentialIntegr_smart(f2, x, y, n) << std::endl;
        std::cout << "seq time: " << clock() - t << std::endl;
    }
}
TEST(Trapezoid_Rule, Test_f1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int dimension = 3;
    double error = 0.001;
    int t;
    std::vector<double> x(dimension);
    std::vector<double> y(dimension);
    x = { 0, 0, 0 };
    y = { 2, 2, 2 };
    const int n = 50;
    t = clock();
    double res = trapezoidParallelRule(f1, x, y, n);
    if (rank == 0) {
        std::cout << "parallel res: " << res << std::endl;
        std::cout << "parallel time: " << clock() - t << std::endl;
        ASSERT_EQ(abs(res - 8.0) < error, 1);
        t = clock();
        std::cout << "seq res: " << SequentialIntegr_smart(f1, x, y, n) << std::endl;
        std::cout << "seq time: " << clock() - t << std::endl;
    }
}

TEST(Trapezoid_Rule, Test_f0) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int dimension = 3;
    double error = 0.001;
    int t;
    std::vector<double> x(dimension);
    std::vector<double> y(dimension);
    x = { 0, 0, 0 };
    y = { 2, 2, 2 };
    const int n = 50;
    t = clock();
    double res = trapezoidParallelRule(f0, x, y, n);
    if (rank == 0) {
        std::cout << "parallel res: " << res << std::endl;
        std::cout << "parallel time: " << clock() - t << std::endl;
        ASSERT_EQ(abs(res - 0.0) < error, 1);
        t = clock();
        std::cout << "seq res: " << SequentialIntegr_smart(f0, x, y, n) << std::endl;
        std::cout << "seq time: " << clock() - t << std::endl;
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
