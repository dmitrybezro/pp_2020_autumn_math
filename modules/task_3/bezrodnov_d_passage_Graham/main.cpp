// Copyright 2020 Bezrodnov Dmitry
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./passage_Graham.h"

TEST(Parallel_Operations_MPI, Cloud_10_point) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    int size = 10;
    std::vector<point> cloud;
    if (RANK == 0) {
    	cloud = getRandomCloud(size);
    }

    std::vector<int> list_parall = ParallelPassageGraham(cloud);

    if (RANK == 0) {
        std::vector<int> list_sequen = SequentialPassageGraham(cloud);
        ASSERT_EQ(list_parall, list_sequen);
    }
}

TEST(Parallel_Operations_MPI, Cloud_57_point) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    int size = 57;
    std::vector<point> cloud;
    if (RANK == 0) {
    	cloud = getRandomCloud(size);
    }

    std::vector<int> list_parall = ParallelPassageGraham(cloud);

    if (RANK == 0) {
        std::vector<int> list_sequen = SequentialPassageGraham(cloud);
        ASSERT_EQ(list_parall, list_sequen);
    }
}

TEST(Parallel_Operations_MPI, Cloud_101_point) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    int size = 101;
    std::vector<point> cloud;
    if (RANK == 0) {
    	cloud = getRandomCloud(size);
    }

    std::vector<int> list_parall = ParallelPassageGraham(cloud);

    if (RANK == 0) {
        std::vector<int> list_sequen = SequentialPassageGraham(cloud);
        ASSERT_EQ(list_parall, list_sequen);
    }
}

TEST(Parallel_Operations_MPI, Cloud_203_point) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    int size = 203;
    std::vector<point> cloud;
    if (RANK == 0) {
    	cloud = getRandomCloud(size);
    }

    std::vector<int> list_parall = ParallelPassageGraham(cloud);

    if (RANK == 0) {
        std::vector<int> list_sequen = SequentialPassageGraham(cloud);
        ASSERT_EQ(list_parall, list_sequen);
    }
}

TEST(Parallel_Operations_MPI, Cloud_256_point) {
    int RANK;
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    int size = 256;
    std::vector<point> cloud;
    if (RANK == 0) {
    	cloud = getRandomCloud(size);
    }

    std::vector<int> list_parall = ParallelPassageGraham(cloud);

    if (RANK == 0) {
        std::vector<int> list_sequen = SequentialPassageGraham(cloud);
        ASSERT_EQ(list_parall, list_sequen);
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