// Copyright 2020 Lebedev Andrew
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <random>
#include <vector>
#include "./my_gather.h"

TEST(my_gather, Test_int) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int count = 10;
    int dest = 1;
    std::vector<int> sbuf(count);
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(NULL)));
    for (int i = 0; i < count; i++) {
        sbuf[i] = gen() % 100;
    }
    double start, stop;
    std::vector<int> rbuf1(count * size), rbuf2(count * size);
    start = MPI_Wtime();
    MPI_Gather(&sbuf[0], count, MPI_INT, &rbuf1[0], count, MPI_INT, dest, MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (rank == dest) {
        printf("MPI_gather time: %f\n", stop - start);
    }
    start = MPI_Wtime();
    my_gather(&sbuf[0], count, MPI_INT, &rbuf2[0], count, MPI_INT, dest, MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (rank == dest) {
        printf("my_gather time: %f\n", stop - start);
        for (int i = 0; i < count * size; i++) {
            EXPECT_EQ(rbuf1[i], rbuf2[i]);
        }
    }
}

TEST(my_gather, Test_0_size) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int count = 0;
    int dest = 1;
    std::vector<int> sbuf(count);
    double start, stop;
    std::vector<int> rbuf1(count * size), rbuf2(count * size);
    start = MPI_Wtime();
    MPI_Gather(&sbuf[0], count, MPI_INT, &rbuf1[0], count, MPI_INT, dest, MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (rank == 0) {
        printf("MPI_gather time: %f\n", stop - start);
    }
    start = MPI_Wtime();
    my_gather(&sbuf[0], count, MPI_INT, &rbuf2[0], count, MPI_INT, dest, MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (rank == dest) {
        printf("my_gather time: %f\n", stop - start);
        for (int i = 0; i < count * size; i++) {
            EXPECT_EQ(rbuf1[i], rbuf2[i]);
        }
    }
}

TEST(my_gather, Test_rand_dest) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(NULL)));
    int count = 10;
    int dest = gen() % size;
    std::vector<int> sbuf(count);
    for (int i = 0; i < count; i++) {
        sbuf[i] = gen() % 100;
    }
    double start, stop;
    std::vector<int> rbuf1(count * size), rbuf2(count * size);
    start = MPI_Wtime();
    MPI_Gather(&sbuf[0], count, MPI_INT, &rbuf1[0], count, MPI_INT, dest, MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (rank == 0) {
        printf("MPI_gather time: %f\n", stop - start);
    }
    start = MPI_Wtime();
    my_gather(&sbuf[0], count, MPI_INT, &rbuf2[0], count, MPI_INT, dest, MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (rank == dest) {
        printf("my_gather time: %f\n", stop - start);
        for (int i = 0; i < count * size; i++) {
            EXPECT_EQ(rbuf1[i], rbuf2[i]);
        }
    }
}

TEST(my_gather, Test_float) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int count = 10;
    int dest = 1;
    std::vector<float> sbuf(count);
    for (int i = 0; i < count; i++) {
        sbuf[i] = i + size * rank / 3;
    }
    double start, stop;
    std::vector<float> rbuf1(count * size), rbuf2(count * size);
    start = MPI_Wtime();
    MPI_Gather(&sbuf[0], count, MPI_FLOAT, &rbuf1[0], count, MPI_FLOAT, dest, MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (rank == 0) {
        printf("MPI_gather time: %f\n", stop - start);
    }
    start = MPI_Wtime();
    my_gather(&sbuf[0], count, MPI_FLOAT, &rbuf2[0], count, MPI_FLOAT, dest, MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (rank == dest) {
        printf("my_gather time: %f\n", stop - start);
        for (int i = 0; i < count * size; i++) {
            EXPECT_EQ(rbuf1[i], rbuf2[i]);
        }
    }
}

TEST(my_gather, Test_double) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int count = 10;
    int dest = 1;
    std::vector<double> sbuf(count);
    for (int i = 0; i < count; i++) {
        sbuf[i] = i + size * rank / 3;
    }
    double start, stop;
    std::vector<double> rbuf1(count * size), rbuf2(count * size);
    start = MPI_Wtime();
    MPI_Gather(&sbuf[0], count, MPI_DOUBLE, &rbuf1[0], count, MPI_DOUBLE, dest, MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (rank == 0) {
        printf("MPI_gather time: %f\n", stop - start);
    }
    start = MPI_Wtime();
    my_gather(&sbuf[0], count, MPI_DOUBLE, &rbuf2[0], count, MPI_DOUBLE, dest, MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (rank == dest) {
        printf("my_gather time: %f\n", stop - start);
        for (int i = 0; i < count * size; i++) {
            EXPECT_EQ(rbuf1[i], rbuf2[i]);
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
