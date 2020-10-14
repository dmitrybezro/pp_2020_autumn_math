// Copyright 2020 Kolesnikov Gleb
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_1/kolesnikov_g_alphabetic_count/aplhabetic_count.h"

std::vector<char> getRandomString(int size) {
    const std::string CHARACTERS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    srand(time(NULL));
    std::vector<char> random_string;
    for (int i = 0; i <size; i++) {
        random_string.push_back(CHARACTERS[rand() % (CHARACTERS.size() - 1)]);
     }
    return random_string;
}

int getSequentialCount(std::vector<char> str) {
    int counter = 0;
    for (int i = 0; i< str.size(); i++) {
        if (str[i] >= 'A' && str[i] <= 'Z' || str[i] >= 'a' && str[i] <= 'z')
            counter++;
    }
    return counter;
}

int getParallelCount(std::vector<char> global_str,
    int vector_size) {
    int size;
    int rank;
    int global_counter = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int delta = vector_size / size;
    if (rank == 0) {
        for (int process = 1; process < size; process++) {
            MPI_Send(&global_str[0] + process * delta, delta, MPI_CHAR, process, 0, MPI_COMM_WORLD);
        }
    }

    std::vector<char> local_str(delta);
    if (rank == 0) {
        local_str = std::vector<char>(global_str.begin(), global_str.begin() + delta);
    }else {
        MPI_Status status;
        MPI_Recv(&local_str[0], delta, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
    }
    int local_counter = getSequentialCount(local_str);
    MPI_Reduce(&local_counter, &global_counter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_counter;
}
