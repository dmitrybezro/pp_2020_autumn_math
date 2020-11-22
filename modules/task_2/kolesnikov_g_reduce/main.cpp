// Copyright 2020 Kolesnikov Gleb
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <random>
#include <vector>
#include <ctime>
#include "./reduce.h"

TEST(MPI_My_Reduce, Int_Sum) {
    int rank;
    int size;
    std::mt19937 gen;
    int count = 10;
    std::vector<int> sendbuf(count);
    std::vector<int> recvbuf(count);
    double t1;
    std::vector<int> check(count);
    gen.seed(static_cast<unsigned int>(time(NULL)));
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int i = 0; i < count; i++) {
        sendbuf[i] = gen() % 13;
    }
    if (rank == 0) {
        t1 = MPI_Wtime();
    }
    MPI_Reduce(&sendbuf[0], &check[0], count, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " origin reduce time" << std::endl;
        t1 = MPI_Wtime();
    }
    MPI_Reduce_My_Own(&sendbuf[0], &recvbuf[0], count, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " my reduce time" << std::endl;
        ASSERT_EQ(recvbuf, check);
    }
}

TEST(MPI_My_Reduce, Int_Min) {
    int rank;
    int size;
    std::mt19937 gen;
    int count = 10;
    std::vector<int> sendbuf(count);
    std::vector<int> recvbuf(count);
    double t1;
    std::vector<int> check(count);
    gen.seed(static_cast<unsigned int>(time(NULL)));
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int i = 0; i < count; i++) {
        sendbuf[i] = gen() % 13;
    }
    if (rank == 0) {
        t1 = MPI_Wtime();
    }
    MPI_Reduce(&sendbuf[0], &check[0], count, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " origin reduce time" << std::endl;
        t1 = MPI_Wtime();
    }
    MPI_Reduce_My_Own(&sendbuf[0], &recvbuf[0], count, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " my reduce time" << std::endl;
        ASSERT_EQ(recvbuf, check);
    }
}

TEST(MPI_My_Reduce, Int_Max) {
    int rank;
    int size;
    std::mt19937 gen;
    int count = 10;
    std::vector<int> sendbuf(count);
    std::vector<int> recvbuf(count);
    double t1;
    std::vector<int> check(count);
    gen.seed(static_cast<unsigned int>(time(NULL)));
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int i = 0; i < count; i++) {
        sendbuf[i] = gen() % 13;
    }
    if (rank == 0) {
        t1 = MPI_Wtime();
    }
    MPI_Reduce(&sendbuf[0], &check[0], count, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " origin reduce time" << std::endl;
        t1 = MPI_Wtime();
    }
    MPI_Reduce_My_Own(&sendbuf[0], &recvbuf[0], count, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " my reduce time" << std::endl;
        ASSERT_EQ(recvbuf, check);
    }
}

TEST(MPI_My_Reduce, Int_Prod) {
    int rank;
    int size;
    std::mt19937 gen;
    int count = 10;
    std::vector<int> sendbuf(count);
    std::vector<int> recvbuf(count);
    double t1;
    std::vector<int> check(count);
    gen.seed(static_cast<unsigned int>(time(NULL)));
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int i = 0; i < count; i++) {
        sendbuf[i] = gen() % 13;
    }
    if (rank == 0) {
        t1 = MPI_Wtime();
    }
    MPI_Reduce(&sendbuf[0], &check[0], count, MPI_INT, MPI_PROD, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " origin reduce time" << std::endl;
        t1 = MPI_Wtime();
    }
    MPI_Reduce_My_Own(&sendbuf[0], &recvbuf[0], count, MPI_INT, MPI_PROD, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " my reduce time" << std::endl;
        ASSERT_EQ(recvbuf, check);
    }
}

TEST(MPI_My_Reduce, Float_Min) {
    int rank;
    int size;
    std::mt19937 gen;
    int count = 3;
    std::vector<float> sendbuf(count);
    std::vector<float> recvbuf(count);
    double t1;
    std::vector<float> check(count);
    gen.seed(static_cast<unsigned int>(time(NULL)));
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int i = 0; i < count; i++) {
        sendbuf[i] = gen();
    }
    if (rank == 0) {
        t1 = MPI_Wtime();
    }
    MPI_Reduce(&sendbuf[0], &check[0], count, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " origin reduce time" << std::endl;
        t1 = MPI_Wtime();
    }
    MPI_Reduce_My_Own(&sendbuf[0], &recvbuf[0], count, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " my reduce time" << std::endl;
        ASSERT_EQ(recvbuf, check);
    }
}
TEST(MPI_My_Reduce, Double_Sum) {
    int rank;
    int size;
    std::mt19937 gen;
    int count = 10;
    std::vector<double> sendbuf(count);
    std::vector<double> recvbuf(count);
    double t1;
    std::vector<double> check(count);
    gen.seed(static_cast<unsigned int>(time(NULL)));
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int i = 0; i < count; i++) {
        sendbuf[i] = gen();
    }
    if (rank == 0) {
        t1 = MPI_Wtime();
    }
    MPI_Reduce(&sendbuf[0], &check[0], count, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " origin reduce time" << std::endl;
        t1 = MPI_Wtime();
    }
    MPI_Reduce_My_Own(&sendbuf[0], &recvbuf[0], count, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " my reduce time" << std::endl;
        ASSERT_EQ(recvbuf, check);
    }
}
TEST(MPI_My_Reduce, Double_Max_with_root) {
    int rank;
    int size;
    std::mt19937 gen;
    int count = 10;
    std::vector<double> sendbuf(count);
    std::vector<double> recvbuf(count);
    double t1;
    int root;
    std::vector<double> check(count);
    gen.seed(static_cast<unsigned int>(time(NULL)));
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1) {
        root = 0;
    } else {
        root = size / 2;
    }
    for (int i = 0; i < count; i++) {
        sendbuf[i] = gen();
    }
    if (rank == 0) {
        t1 = MPI_Wtime();
    }
    MPI_Reduce(&sendbuf[0], &check[0], count, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " origin reduce time" << std::endl;
        t1 = MPI_Wtime();
    }
    MPI_Reduce_My_Own(&sendbuf[0], &recvbuf[0], count, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << MPI_Wtime() - t1 << " my reduce time" << std::endl;
    }
    if (rank == root) {
        ASSERT_EQ(recvbuf, check);
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
