// Copyright 2020 Mishin Ilya
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include "./odd_even_shell_sort.h"

TEST(Parallel_Operations_MPI, Test_task_1) {
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> vec(1), vec1(1), vec2(1);
    int N = 512;
    if (rank == 0) {
        vec = getRandomVector(N);
        vec1 = vec;
        vec2 = vec;
    }

    double t1;
    if (rank == 0) {
        t1 = MPI_Wtime();
    }

    parallel_shell_sort(&vec1[0], N);

    if (rank == 0) {
        std::cout << "My Sort Time: " << MPI_Wtime() - t1 << std::endl;
        t1 = MPI_Wtime();
        shell_sort(vec2.begin(), vec2.end(), [](int a, int b) {
            return a < b;
        });
        std::cout << "Sort Time: " << MPI_Wtime() - t1 << std::endl;
        ASSERT_EQ(vec1, vec2);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

TEST(Parallel_Operations_MPI, Test_task_2) {
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> vec(1), vec1(1), vec2(1);
    int N = 1024;
    if (rank == 0) {
        vec = getRandomVector(N);
        vec1 = vec;
        vec2 = vec;
    }

    double t1;
    if (rank == 0) {
        t1 = MPI_Wtime();
    }

    parallel_shell_sort(&vec1[0], N);

    if (rank == 0) {
        std::cout << "My Sort Time: " << MPI_Wtime() - t1 << std::endl;
        t1 = MPI_Wtime();
        shell_sort(vec2.begin(), vec2.end(), [](int a, int b) {
            return a < b;
        });
        std::cout << "Sort Time: " << MPI_Wtime() - t1 << std::endl;
        ASSERT_EQ(vec1, vec2);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

TEST(Parallel_Operations_MPI, Test_task_3) {
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> vec(1), vec1(1), vec2(1);
    int N = 2048;
    if (rank == 0) {
        vec = getRandomVector(N);
        vec1 = vec;
        vec2 = vec;
    }

    double t1;
    if (rank == 0) {
        t1 = MPI_Wtime();
    }

    parallel_shell_sort(&vec1[0], N);

    if (rank == 0) {
        std::cout << "My Sort Time: " << MPI_Wtime() - t1 << std::endl;
        t1 = MPI_Wtime();
        shell_sort(vec2.begin(), vec2.end(), [](int a, int b) {
            return a < b;
        });
        std::cout << "Sort Time: " << MPI_Wtime() - t1 << std::endl;
        ASSERT_EQ(vec1, vec2);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

TEST(Parallel_Operations_MPI, Test_task_4) {
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> vec(1), vec1(1), vec2(1);
    int N = 4096;
    if (rank == 0) {
        vec = getRandomVector(N);
        vec1 = vec;
        vec2 = vec;
    }

    double t1;
    if (rank == 0) {
        t1 = MPI_Wtime();
    }

    parallel_shell_sort(&vec1[0], N);

    if (rank == 0) {
        std::cout << "My Sort Time: " << MPI_Wtime() - t1 << std::endl;
        t1 = MPI_Wtime();
        shell_sort(vec2.begin(), vec2.end(), [](int a, int b) {
            return a < b;
        });
        std::cout << "Sort Time: " << MPI_Wtime() - t1 << std::endl;
        ASSERT_EQ(vec1, vec2);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

TEST(Parallel_Operations_MPI, Test_task_5) {
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> vec(1), vec1(1), vec2(1);
    int N = 8192;
    if (rank == 0) {
        vec = getRandomVector(N);
        vec1 = vec;
        vec2 = vec;
    }

    double t1;
    if (rank == 0) {
        t1 = MPI_Wtime();
    }

    parallel_shell_sort(&vec1[0], N);

    if (rank == 0) {
        std::cout << "My Sort Time: " << MPI_Wtime() - t1 << std::endl;
        t1 = MPI_Wtime();
        shell_sort(vec2.begin(), vec2.end(), [](int a, int b) {
            return a < b;
        });
        std::cout << "Sort Time: " << MPI_Wtime() - t1 << std::endl;
        ASSERT_EQ(vec1, vec2);
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
