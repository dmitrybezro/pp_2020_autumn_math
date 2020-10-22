// Copyright 2018 Nesterov Alexander
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./matrix_sum.h"

TEST(Parallel_Operations_MPI, Test_Sum) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_mat;
    const int count_row = 5;
	const int count_str = 5;

    if (rank == 0) {
        global_mat = getRandomMat(count_row,count_str);
    }

    std::vector<int> global_sum = getParallelOperations(global_mat, count_row,count_str);

    if (rank == 0) {
        std::vector<int> reference_sum = getSequentialOperations(global_mat,count_row,count_str);
        ASSERT_EQ(reference_sum, global_sum);
    }
}
TEST(Parallel_Operations_MPI, Test_Sum_2) {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::vector<int> global_mat;
	const int count_row = 100;
	const int count_str = 100;

	if (rank == 0) {
		global_mat = getRandomMat(count_row, count_str);
	}

	std::vector<int> global_sum = getParallelOperations(global_mat, count_row, count_str);

	if (rank == 0) {
		std::vector<int> reference_sum = getSequentialOperations(global_mat, count_row, count_str);
		ASSERT_EQ(reference_sum, global_sum);
	}
}
TEST(Parallel_Operations_MPI, Test_Sum_3) {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::vector<int> global_mat;
	const int count_row = 105;
	const int count_str = 213;

	if (rank == 0) {
		global_mat = getRandomMat(count_row, count_str);
	}

	std::vector<int> global_sum = getParallelOperations(global_mat, count_row, count_str);

	if (rank == 0) {
		std::vector<int> reference_sum = getSequentialOperations(global_mat, count_row, count_str);
		ASSERT_EQ(reference_sum, global_sum);
	}
}
TEST(Parallel_Operations_MPI, Test_Sum_4) {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::vector<int> global_mat;
	const int count_row = 1002;
	const int count_str = 2001;

	if (rank == 0) {
		global_mat = getRandomMat(count_row, count_str);
	}

	std::vector<int> global_sum = getParallelOperations(global_mat, count_row, count_str);

	if (rank == 0) {
		std::vector<int> reference_sum = getSequentialOperations(global_mat, count_row, count_str);
		ASSERT_EQ(reference_sum, global_sum);
	}
}
TEST(Parallel_Operations_MPI, Test_Sum_5) {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::vector<int> global_mat;
	const int count_row = 51;
	const int count_str = 43;

	if (rank == 0) {
		global_mat = getRandomMat(count_row, count_str);
	}

	std::vector<int> global_sum = getParallelOperations(global_mat, count_row, count_str);

	if (rank == 0) {
		std::vector<int> reference_sum = getSequentialOperations(global_mat, count_row, count_str);
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
