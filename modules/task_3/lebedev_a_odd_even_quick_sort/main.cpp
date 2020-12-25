// Copyright 2020 Lebedev Andrew

#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <ctime>
#include <iostream>
#include <cmath>
#include "../../../modules/task_3/lebedev_a_odd_even_quick_sort/odd_even_quick_sort.h"

TEST(odd_even_qsort, test_4) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 4;
    double start, finish;
    std::vector<int> vec, vec1;
    if (rank == 0) {
        vec = vec1 = getRandomVector(size);
        start = MPI_Wtime();
    }
    parallel_qsort(vec.data(), size);
    if (rank == 0) {
        finish = MPI_Wtime();
        printf("My time: %lf\n", finish - start);
        start = MPI_Wtime();
        qsort(vec1.data(), size, sizeof(int), compare);
        finish = MPI_Wtime();
        printf("qsort time: %lf\n", finish - start);
        EXPECT_EQ(vec, vec1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

TEST(odd_even_qsort, test_32) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 32;
    double start, finish;
    std::vector<int> vec, vec1;
    if (rank == 0) {
        vec = vec1 = getRandomVector(size);
        start = MPI_Wtime();
    }
    parallel_qsort(vec.data(), size);
    if (rank == 0) {
        finish = MPI_Wtime();
        printf("My time: %lf\n", finish - start);
        start = MPI_Wtime();
        qsort(vec1.data(), size, sizeof(int), compare);
        finish = MPI_Wtime();
        printf("qsort time: %lf\n", finish - start);
        EXPECT_EQ(vec, vec1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

TEST(odd_even_qsort, test_64) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 64;
    double start, finish;
    std::vector<int> vec, vec1;
    if (rank == 0) {
        vec = vec1 = getRandomVector(size);
        start = MPI_Wtime();
    }
    parallel_qsort(vec.data(), size);
    if (rank == 0) {
        finish = MPI_Wtime();
        printf("My time: %lf\n", finish - start);
        start = MPI_Wtime();
        qsort(vec1.data(), size, sizeof(int), compare);
        finish = MPI_Wtime();
        printf("qsort time: %lf\n", finish - start);
        EXPECT_EQ(vec, vec1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

TEST(odd_even_qsort, test_128) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 128;
    double start, finish;
    std::vector<int> vec, vec1;
    if (rank == 0) {
        vec = vec1 = getRandomVector(size);
        start = MPI_Wtime();
    }
    parallel_qsort(vec.data(), size);
    if (rank == 0) {
        finish = MPI_Wtime();
        printf("My time: %lf\n", finish - start);
        start = MPI_Wtime();
        qsort(vec1.data(), size, sizeof(int), compare);
        finish = MPI_Wtime();
        printf("qsort time: %lf\n", finish - start);
        EXPECT_EQ(vec, vec1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

TEST(odd_even_qsort, test_256) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 256;
    double start, finish;
    std::vector<int> vec, vec1;
    if (rank == 0) {
        vec = vec1 = getRandomVector(size);
        start = MPI_Wtime();
    }
    parallel_qsort(vec.data(), size);
    if (rank == 0) {
        finish = MPI_Wtime();
        printf("My time: %lf\n", finish - start);
        start = MPI_Wtime();
        qsort(vec1.data(), size, sizeof(int), compare);
        finish = MPI_Wtime();
        printf("qsort time: %lf\n", finish - start);
        EXPECT_EQ(vec, vec1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
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
