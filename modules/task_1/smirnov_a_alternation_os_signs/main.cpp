// Copyright 2020 Smirnov Alexander
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include "./alternation_os_signs.h"

// #define time101 true

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
TEST(Parallel_Operations_MPI, Test_10000000) {
    int rank;
    double t1, t2;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> global_vec;
    const int count_size_vector = 100000000;
    int global_sum = 0;
    int reference_sum = 0;
    int n = 101;  // count iteration (101)
    std::vector<double> timePar(n, 0);
    std::vector<double> timeSeq(n, 0);

    if (rank == 0) {
        global_vec = getRandomVector(count_size_vector);
    }
#ifdef time101  // PAR
    for (int i = 0; i < n + 1; i++) {
        if (rank == 0)
            t1 = MPI_Wtime();
#endif  // time101 // PAR
        global_sum = getParallelOperations(global_vec, count_size_vector);
#ifdef time101  // PAR
        if (rank == 0) {
            t2 = MPI_Wtime();
            if (i != 0)
                timePar[i - 1] = t2 - t1;
        }
    }
    if (rank == 0) {
        std::sort(timePar.begin(), timePar.end());
        printf("time PAR= %3.20f\n", timePar[timePar.size()/2]);
    }
#endif  // time101 // PAR

    if (rank == 0) {
#ifdef time101  // SEQ
        for (int i = 0; i < n + 1; i++) {
            if (rank == 0)
                t1 = MPI_Wtime();
#endif  // time101 // SEQ
            reference_sum = getSequentialOperations(global_vec);
#ifdef time101  // SEQ
            if (rank == 0) {
                t2 = MPI_Wtime();
                if (i != 0)
                    timePar[i - 1] = t2 - t1;
            }
        }
        if (rank == 0) {
            std::sort(timePar.begin(), timePar.end());
            printf("time SEQ= %3.20f\n", timePar[n / 2]);
        }
#endif  // time101 // SEQ
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
