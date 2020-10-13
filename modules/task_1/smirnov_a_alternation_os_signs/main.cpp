// Copyright 2020 Smirnov Alexander
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./alternation_os_signs.h"

TEST(Parallel_Operations_MPI, Test_0) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_vec;
    const int count_size_vector = 0;

    if (rank == 0) {
        global_vec = getRandomVector(count_size_vector);
    }

    int global_sum = getParallelOperations(global_vec, count_size_vector);

    if (rank == 0) {
        int reference_sum = getSequentialOperations(global_vec);
        ASSERT_EQ(reference_sum, global_sum);
    }
}
TEST(Parallel_Operations_MPI, Test_1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_vec;
    const int count_size_vector = 1;

    if (rank == 0) {
        global_vec = getRandomVector(count_size_vector);
    }

    int global_sum = getParallelOperations(global_vec, count_size_vector);

    if (rank == 0) {
        int reference_sum = getSequentialOperations(global_vec);
        ASSERT_EQ(reference_sum, global_sum);
    }
}
TEST(Parallel_Operations_MPI, Test_10) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_vec;
    const int count_size_vector = 10;

    if (rank == 0) {
        global_vec = getRandomVector(count_size_vector);
    }

    int global_sum = getParallelOperations(global_vec, count_size_vector);

    if (rank == 0) {
        int reference_sum = getSequentialOperations(global_vec);
        ASSERT_EQ(reference_sum, global_sum);
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

    int global_sum = getParallelOperations(global_vec, count_size_vector);

    if (rank == 0) {
        int reference_sum = getSequentialOperations(global_vec);
        ASSERT_EQ(reference_sum, global_sum);
    }
}
TEST(Parallel_Operations_MPI, Test_5000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_vec;
    const int count_size_vector = 5000;

    if (rank == 0) {
        global_vec = getRandomVector(count_size_vector);
    }

    int global_sum = getParallelOperations(global_vec, count_size_vector);

    if (rank == 0) {
        int reference_sum = getSequentialOperations(global_vec);
        ASSERT_EQ(reference_sum, global_sum);
    }
}
TEST(Parallel_Operations_MPI, Test_1000000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_vec;
    const int count_size_vector = 1000000;

    if (rank == 0) {
        global_vec = getRandomVector(count_size_vector);
    }

    int global_sum = getParallelOperations(global_vec, count_size_vector);

    if (rank == 0) {
        int reference_sum = getSequentialOperations(global_vec);
        ASSERT_EQ(reference_sum, global_sum);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    int ProcNum, ProcRank;
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

#ifdef detail
    if (ProcRank == 0) std::cout << ProcNum << std::endl;
#endif
    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
