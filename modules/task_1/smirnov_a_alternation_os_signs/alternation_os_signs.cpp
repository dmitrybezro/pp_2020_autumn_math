// Copyright 2018 Nesterov Alexander
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_1/smirnov_a_alternation_os_signs/alternation_os_signs.h"


std::vector<int> getRandomVector(int sz) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> vec(sz);
    for (int i = 0; i < sz; i++) { vec[i] = gen() % 100 - 50; }
    return vec;
}

int getSequentialOperations(std::vector<int> vec) {
    const int  sz = vec.size();
    if (sz == 0) return 0;
    double t1, t2;
    int change_of_sings = 0;
    t1 = MPI_Wtime();
    for (int i = 0; i < sz - 1; i++) {
        if (vec[i] * vec[i + 1] < 0) change_of_sings++;
    }
    t2 = MPI_Wtime();
    printf("time = %3.20f\n", t2 - t1);
    return change_of_sings;
}

int getParallelOperations(std::vector<int> global_vec,
    int count_size_vector) {
    double t1, t2;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int change_of_sings = 0;

    if (count_size_vector == 0) return 0;
    if (rank == 0) {
        t1 = MPI_Wtime();
        for (int proc = 1; proc < size; proc++) {
            MPI_Send(&global_vec[0], count_size_vector,
                MPI_INT, proc, 0, MPI_COMM_WORLD);
        }
    }

    std::vector<int> local_vec(count_size_vector);
    if (rank == 0) {
        local_vec = std::vector<int>(global_vec.begin(),
            global_vec.begin() + count_size_vector);
    } else {
        MPI_Status status;
        MPI_Recv(&local_vec[0], count_size_vector, MPI_INT, ,
            0, 0, MPI_COMM_WORLD, &status);
    }

    int local_count = 0;
    for (int i = rank; i < (count_size_vector - 1); i += size) {
        bool positive;
        if (local_vec[i] * local_vec[i + 1] < 0) local_count++;
    }
    MPI_Reduce(&local_count, &change_of_sings, 1, MPI_INT,
        MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        t2 = MPI_Wtime();
        printf("time = %3.20f\n", t2 - t1);
    }
    return change_of_sings;
}
