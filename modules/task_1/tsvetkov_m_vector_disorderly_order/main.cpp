// Copyright 2020 Tsvetkov Maxim
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./vector_disorderly_order.h"



TEST(Parallel_Operations_MPI, Test_55) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec;
    const int count_size = 55;

    if (rank == 0) {
        vec = getRandomVector(count_size);
    }

    int global_sum = get_parallel_operations(vec, count_size);

    if (rank == 0) {
        int reference_sum = get_sequential_operations(vec);
        ASSERT_EQ(reference_sum, global_sum);
    }
}

TEST(Parallel_Operations_MPI, Test_15) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec;
    const int count_size = 15;

    if (rank == 0) {
        vec = getRandomVector(count_size);
    }

    int global_sum = get_parallel_operations(vec, count_size);

    if (rank == 0) {
        int reference_sum = get_sequential_operations(vec);
        ASSERT_EQ(reference_sum, global_sum);
    }
}

TEST(Parallel_Operations_MPI, Test_23) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec;
    const int count_size = 23;

    if (rank == 0) {
        vec = getRandomVector(count_size);
    }

    int global_sum = get_parallel_operations(vec, count_size);

    if (rank == 0) {
        int reference_sum = get_sequential_operations(vec);
        ASSERT_EQ(reference_sum, global_sum);
    }
}

TEST(Parallel_Operations_MPI, Test_254) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec;
    const int count_size = 254;

    if (rank == 0) {
        vec = getRandomVector(count_size);
    }

    int global_sum = get_parallel_operations(vec, count_size);

    if (rank == 0) {
        int reference_sum = get_sequential_operations(vec);
        ASSERT_EQ(reference_sum, global_sum);
    }
}


TEST(Parallel_Operations_MPI, Test_512) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec;
    const int count_size = 512;

    if (rank == 0) {
        vec = getRandomVector(count_size);
    }

    int global_sum = get_parallel_operations(vec, count_size);

    if (rank == 0) {
        int reference_sum = get_sequential_operations(vec);
        ASSERT_EQ(reference_sum, global_sum);
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
