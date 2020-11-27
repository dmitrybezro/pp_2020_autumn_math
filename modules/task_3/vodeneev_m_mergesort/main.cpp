// Copyright 2020 Vodeneev Mikhail
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./ops_mpi.h"

TEST(Parallel_Operations_MPI, Test_1) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec(10), vec_res(10);
    int size_vec = vec.size();
    if (rank == 0) {
        int k = vec.size();
        for (int i = 0; i < size_vec; i++) {
            vec[i] = k;
            k--;
        }
    }
    double time1 = MPI_Wtime();
    vec_res = mergesort(vec);
    double time2 = MPI_Wtime();
    double time_res1 = time2 - time1;
    if (rank == 0) {
        double time3 = MPI_Wtime();
        qsort(vec.data(), vec.size(), sizeof(vec[0]), compare);
        double time4 = MPI_Wtime();
        double time_res2 = time4 - time3;
        std::cout << time_res1 << ' ' << time_res2 << std :: endl;
        ASSERT_EQ(vec_res[3], vec[3]);
    }
}

TEST(Parallel_Operations_MPI, Test_2) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec(100), vec_res(100);
    int size_vec = vec.size();
    if (rank == 0) {
        int k = vec.size();
        for (int i = 0; i < size_vec; i++) {
            vec[i] = k;
            k--;
        }
    }
    double time1 = MPI_Wtime();
    vec_res = mergesort(vec);
    double time2 = MPI_Wtime();
    double time_res1 = time2 - time1;
    if (rank == 0) {
        double time3 = MPI_Wtime();
        qsort(vec.data(), vec.size(), sizeof(vec[0]), compare);
        double time4 = MPI_Wtime();
        double time_res2 = time4 - time3;
        std::cout << time_res1 << ' ' << time_res2 << std::endl;
        ASSERT_EQ(vec_res[30], vec[30]);
    }
}

TEST(Parallel_Operations_MPI, Test_3) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec(10000), vec_res(10000);
    int size_vec = vec.size();
    if (rank == 0) {
        int k = vec.size();
        for (int i = 0; i < size_vec; i++) {
            vec[i] = k;
            k--;
        }
    }
    double time1 = MPI_Wtime();
    vec_res = mergesort(vec);
    double time2 = MPI_Wtime();
    double time_res1 = time2 - time1;
    if (rank == 0) {
        double time3 = MPI_Wtime();
        qsort(vec.data(), vec.size(), sizeof(vec[0]), compare);
        double time4 = MPI_Wtime();
        double time_res2 = time4 - time3;
        std::cout << time_res1 << ' ' << time_res2 << std::endl;
        ASSERT_EQ(vec_res[30], vec[30]);
    }
}

TEST(Parallel_Operations_MPI, Test_4) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec(100000), vec_res(100000);
    int size_vec = vec.size();
    if (rank == 0) {
        int k = vec.size();
        for (int i = 0; i < size_vec; i++) {
            vec[i] = k;
            k--;
        }
    }
    double time1 = MPI_Wtime();
    vec_res = mergesort(vec);
    double time2 = MPI_Wtime();
    double time_res1 = time2 - time1;
    if (rank == 0) {
        double time3 = MPI_Wtime();
        qsort(vec.data(), vec.size(), sizeof(vec[0]), compare);
        double time4 = MPI_Wtime();
        double time_res2 = time4 - time3;
        std::cout << time_res1 << ' ' << time_res2 << std::endl;
        ASSERT_EQ(vec_res[400], vec[400]);
    }
}

TEST(Parallel_Operations_MPI, Test_5) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec(10000000), vec_res(10000000);
    int size_vec = vec.size();
    if (rank == 0) {
        int k = vec.size();
        for (int i = 0; i < size_vec; i++) {
            vec[i] = k;
            k--;
        }
    }
    double time1 = MPI_Wtime();
    vec_res = mergesort(vec);
    double time2 = MPI_Wtime();
    double time_res1 = time2 - time1;
    if (rank == 0) {
        double time3 = MPI_Wtime();
        qsort(vec.data(), vec.size(), sizeof(vec[0]), compare);
        double time4 = MPI_Wtime();
        double time_res2 = time4 - time3;
        std::cout << time_res1 << ' ' << time_res2 << std::endl;
        ASSERT_EQ(vec_res[4000], vec[4000]);
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
