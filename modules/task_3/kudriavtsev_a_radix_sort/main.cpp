// Copyright 2020 Kudriavtsev Alexander
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./ops_mpi.h"

TEST(Parallel_Operations_MPI, Test_Seq_20) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec1;
    const int n = 20;
    if (rank == 0) {
        std::vector<int> vec2 = vec1 = getRandomVector(n);
        double t_b = MPI_Wtime();
        radixSort(vec1.data(), n);
        double t_e = MPI_Wtime();
        std::cout << "Sequential time: " << t_e - t_b << std::endl;
        std::sort(vec2.begin(), vec2.end());
        for (int i = 0; i < n; ++i)
            ASSERT_EQ(vec2[i], vec1[i]);
    }
}

TEST(Parallel_Operations_MPI, Test_Seq_2000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec1;
    const int n = 20;
    if (rank == 0) {
        std::vector<int> vec2 = vec1 = getRandomVector(n);
        double t_b = MPI_Wtime();
        radixSort(vec1.data(), n);
        double t_e = MPI_Wtime();
        std::cout << "Sequential time: " << t_e - t_b << std::endl;
        std::sort(vec2.begin(), vec2.end());
        for (int i = 0; i < n; ++i)
            ASSERT_EQ(vec2[i], vec1[i]);
    }
}

TEST(Parallel_Operations_MPI, Test_Seq_20000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec1;
    const int n = 20;
    if (rank == 0) {
        std::vector<int> vec2 = vec1 = getRandomVector(n);
        double t_b = MPI_Wtime();
        radixSort(vec1.data(), n);
        double t_e = MPI_Wtime();
        std::cout << "Sequential time: " << t_e - t_b << std::endl;
        std::sort(vec2.begin(), vec2.end());
        for (int i = 0; i < n; ++i)
            ASSERT_EQ(vec2[i], vec1[i]);
    }
}

TEST(Parallel_Operations_MPI, Test_Par_20) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec1, vec2;
    const int n = 20;
    if (rank == 0) {
        vec2 = vec1 = getRandomVector(n);
    }
    parallelRadixSort(vec1.data(), n);
    if (rank == 0) {
        std::sort(vec2.begin(), vec2.end());
        for (int i = 0; i < n; ++i)
            ASSERT_EQ(vec2[i], vec1[i]);
    }
}

TEST(Parallel_Operations_MPI, Test_Par_2000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec1, vec2;
    const int n = 2000;
    if (rank == 0) {
        vec2 = vec1 = getRandomVector(n);
    }
    parallelRadixSort(vec1.data(), n);
    if (rank == 0) {
        std::sort(vec2.begin(), vec2.end());
        for (int i = 0; i < n; ++i)
            ASSERT_EQ(vec2[i], vec1[i]);
    }
}

TEST(Parallel_Operations_MPI, Test_Par_3041) {  // Simple Number
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec1, vec2;
    const int n = 3041;
    if (rank == 0) {
        vec2 = vec1 = getRandomVector(n);
    }
    parallelRadixSort(vec1.data(), n);
    if (rank == 0) {
        std::sort(vec2.begin(), vec2.end());
        for (int i = 0; i < n; ++i)
            ASSERT_EQ(vec2[i], vec1[i]);
    }
}

TEST(Parallel_Operations_MPI, Test_Par_43201) {  // Big Simple Number
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec1, vec2;
    const int n = 43201;
    if (rank == 0) {
        vec2 = vec1 = getRandomVector(n);
    }
    parallelRadixSort(vec1.data(), n);
    if (rank == 0) {
        std::sort(vec2.begin(), vec2.end());
        for (int i = 0; i < n; ++i)
            ASSERT_EQ(vec2[i], vec1[i]);
    }
}

TEST(Parallel_Operations_MPI, Test_Par_1000000) {  // Big Number 1
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec1, vec2;
    const int n = 1000000;
    if (rank == 0) {
        vec2 = vec1 = getRandomVector(n);
    }
    parallelRadixSort(vec1.data(), n);
    if (rank == 0) {
        double t_b = MPI_Wtime();
        radixSort(vec2.data(), n);
        double t_e = MPI_Wtime();
        std::cout << "Sequential time: " << t_e - t_b << std::endl;
        for (int i = 0; i < n; ++i)
            ASSERT_EQ(vec2[i], vec1[i]);
    }
}

TEST(Parallel_Operations_MPI, Test_Par_5000000) {  // Big Number 2
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec1, vec2;
    const int n = 5000000;
    if (rank == 0) {
        vec2 = vec1 = getRandomVector(n);
    }
    parallelRadixSort(vec1.data(), n);
    if (rank == 0) {
        double t_b = MPI_Wtime();
        radixSort(vec2.data(), n);
        double t_e = MPI_Wtime();
        std::cout << "Sequential time: " << t_e - t_b << std::endl;
        for (int i = 0; i < n; ++i)
            ASSERT_EQ(vec2[i], vec1[i]);
    }
}

TEST(Parallel_Operations_MPI, Test_Par_10000000) {  // Big Number 3
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec1, vec2;
    const int n = 10000000;
    if (rank == 0) {
        vec2 = vec1 = getRandomVector(n);
    }
    parallelRadixSort(vec1.data(), n);
    if (rank == 0) {
        double t_b = MPI_Wtime();
        radixSort(vec2.data(), n);
        double t_e = MPI_Wtime();
        std::cout << "Sequential time: " << t_e - t_b << std::endl;
        for (int i = 0; i < n; ++i)
            ASSERT_EQ(vec2[i], vec1[i]);
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
