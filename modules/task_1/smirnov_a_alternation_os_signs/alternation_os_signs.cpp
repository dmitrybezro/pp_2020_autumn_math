    // Copyright 2020 Smirnov Aleksandr
    #include <mpi.h>
    #include <vector>
    #include <string>
    #include <random>
    #include <ctime>
    #include <algorithm>
    #include "./alternation_os_signs.h"

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
    printf("time seq= %3.20f\n", t2 - t1);
    return change_of_sings;
    }

    int getParallelOperations(std::vector<int> global_vec,
    int count_size_vector) {
    double t1, t2;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int change_of_sings = 0;

    int* local_vec = new int[count_size_vector];
    if (count_size_vector == 0) return 0;
    if (rank == 0) {
        for (int i = 0; i < count_size_vector; i++)
            local_vec[i] = global_vec[i];
        t1 = MPI_Wtime();
    }
    MPI_Bcast(local_vec, count_size_vector, MPI_INT, 0, MPI_COMM_WORLD);

    int local_count = 0;
    for (int i = rank; i < (count_size_vector - 1); i += size) {
        if (local_vec[i] * local_vec[i + 1] < 0) local_count++;
    }
    MPI_Reduce(&local_count, &change_of_sings, 1, MPI_INT,
        MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        t2 = MPI_Wtime();
        printf("time par= %3.20f\n", t2 - t1);
    }
    return change_of_sings;
    }
