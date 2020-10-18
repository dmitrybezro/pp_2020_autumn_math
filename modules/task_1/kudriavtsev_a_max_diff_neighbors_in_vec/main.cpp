// Copyright 2020 Kudriavtsev Alexander
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./ops_mpi.h"

TEST(Parallel_Operations_MPI, Test_20) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_vec;
    const int count_size_vector = 20;
    if (rank == 0) {
        global_vec = getRandomVector(count_size_vector);
    }
    int global_res = getParallelOperations(global_vec, count_size_vector);
    if (rank == 0) {
        int reference_res = getSequentialOperations(global_vec, count_size_vector);
        ASSERT_EQ(reference_res, global_res);
    }
}

TEST(Parallel_Operations_MPI, Test_100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_vec;
    const int count_size_vector = 100;
    if (rank == 0) {
        global_vec = getRandomVector(count_size_vector);
    }
    int global_res = getParallelOperations(global_vec, count_size_vector);
    if (rank == 0) {
        int reference_res = getSequentialOperations(global_vec, count_size_vector);
        ASSERT_EQ(reference_res, global_res);
    }
}

TEST(Parallel_Operations_MPI, Test_1000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_vec;
    const int count_size_vector = 1000;
    if (rank == 0) {
        global_vec = getRandomVector(count_size_vector);
    }
    int global_res = getParallelOperations(global_vec, count_size_vector);
    if (rank == 0) {
        int reference_res = getSequentialOperations(global_vec, count_size_vector);
        ASSERT_EQ(reference_res, global_res);
    }
}

TEST(Parallel_Operations_MPI, Test_1001) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_vec;
    const int count_size_vector = 1001;
    if (rank == 0) {
        global_vec = getRandomVector(count_size_vector);
    }
    int global_res = getParallelOperations(global_vec, count_size_vector);
    if (rank == 0) {
        int reference_res = getSequentialOperations(global_vec, count_size_vector);
        std::cout << global_res << " " << reference_res << std::endl;
        ASSERT_EQ(reference_res, global_res);
    }
}

TEST(Parallel_Operations_MPI, Test_2957) {  // Simple Number
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_vec;
    const int count_size_vector = 2957;
    if (rank == 0) {
        global_vec = getRandomVector(count_size_vector);
    }
    int global_res = getParallelOperations(global_vec, count_size_vector);
    if (rank == 0) {
        int reference_res = getSequentialOperations(global_vec, count_size_vector);
        std::cout << global_res << " " << reference_res << std::endl;
        ASSERT_EQ(reference_res, global_res);
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
