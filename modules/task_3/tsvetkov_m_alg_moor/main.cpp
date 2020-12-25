// Copyright 2020 Tsvetkov Maxim
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>

#include "./moors_algorithm.h"

TEST(Moors_Algorithm_MPI, Cant_Find_If_Incorrect_Source_Seq) {
    int rank;
    int n = 3;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> g(n * n);
    if (rank == 0)
        g = getRandomGraph(n);
    std::vector<int> ans(n);
    if (rank == 0) {
        ASSERT_ANY_THROW(Moor_alg(g, (-1)));
    }
}

TEST(Moors_Algorithm_MPI, Can_Find_If_Correct_Source_Seq) {
    int rank;
    int n = 3;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> g(n * n);
    if (rank == 0)
        g = getRandomGraph(n);
    std::vector<int> ans(n);
    if (rank == 0) {
        ASSERT_NO_THROW(Moor_alg(g, (1)));
    }
}

TEST(Moors_Algorithm_MPI, Test_On_Size_55) {
    int rank;
    int n = 55;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> g(n * n);
    if (rank == 0)
        g = getRandomGraph(n);
    std::vector<int> ans1(n);
    std::vector<int> ans2(n);
    ans1 = ParallelMoor(g, 0);
    if (rank == 0) {
        ans2 = Moor_alg(g, 0);
        for (int i = 0; i < n; ++i) {
            ASSERT_EQ(ans2[i], ans1[i]);
        }
    }
}

TEST(Moors_Algorithm_MPI, Test_On_Size_110) {
    int rank;
    int n = 110;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> g(n * n);
    if (rank == 0)
        g = getRandomGraph(n);
    std::vector<int> ans1(n);
    std::vector<int> ans2(n);
    ans1 = ParallelMoor(g, 0);
    if (rank == 0) {
        ans2 = Moor_alg(g, 0);
        for (int i = 0; i < n; ++i) {
            ASSERT_EQ(ans2[i], ans1[i]);
        }
    }
}

TEST(Moors_Algorithm_MPI, Test_On_Size_5) {
    int rank;
    int n = 5;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> g(n * n);
    if (rank == 0)
        g = getRandomGraph(n);
    std::vector<int> ans1(n);
    std::vector<int> ans2(n);
    ans1 = ParallelMoor(g, 0);
    if (rank == 0) {
        ans2 = Moor_alg(g, 0);
        for (int i = 0; i < n; ++i) {
            ASSERT_EQ(ans2[i], ans1[i]);
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
