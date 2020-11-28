// Copyright 2020 Mishin Ilya
#include <mpi.h>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>
#include <cstring>
#include "../../../modules/task_3/mishin_i_odd_even_shell_sort/odd_even_shell_sort.h"

template <class Item>
void batchersort(Item* a, int l, int r) {
int N = r-l;
        for (int p = 1; p < N; p *= 2)
            for (int k = p; k > 0; k /= 2)
                for (int j = k % p; j + k < N; j += (k + k))
                    for (int i = 0; i < N - j - k; i++)
                        if ((j + i) / (p + p) == (j + i + k) / (p + p))
                            if (a[l + j + i] > a[l + j + i + k])
                                std::swap(a[l + j + i], a[l + j + i + k]);
}

int* Batcher_merge(int* array1, int* array2, int size1, int size2) {
    int* result = new int[size1 + size2];
    std::memcpy(result, array1, size1 * sizeof(int));
    std::memcpy(result + size1, array2, size2 * sizeof(int));

    batchersort(result, 0,  size1 + size2);
    delete[] array1;
    delete[] array2;
    return result;
}

int* parallel_shell_sort(int* container, int array_size) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n = array_size;
    int current_array_size;
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank < n % size)
        current_array_size = n / size + 1;
    else
        current_array_size = n / size;

    int* sendcounts = new int[size];
    int* displs = new int[size];

    if (rank == 0) {
        int dif = 0;
        for (int i = 0; i < size; i++) {
            sendcounts[i] = (i < n% size) ? (n / size + 1) : (n / size);
            displs[i] = dif;
            dif += (i < n% size) ? (n / size + 1) : (n / size);
        }
    }

    int* x = new int[current_array_size];

    MPI_Scatterv(container, sendcounts, displs, MPI_INT, x, current_array_size, MPI_INT, 0, MPI_COMM_WORLD);

    shell_sort(x, x + current_array_size, [](int a, int b) {
        return a < b;
    });

    int s = size, m = 1;

    while (s > 1) {
        s = s / 2 + s % 2;
        if ((rank - m) % (2 * m) == 0) {
            MPI_Send(&current_array_size, 1, MPI_INT, rank - m, 0, MPI_COMM_WORLD);
            MPI_Send(x, current_array_size, MPI_INT, rank - m, 0, MPI_COMM_WORLD);
        }
        if ((rank % (2 * m) == 0) && (size - rank > m)) {
            MPI_Status status;
            int input_array_size;
            MPI_Recv(&input_array_size, 1, MPI_INT, rank + m, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int* y = new int[input_array_size];
            MPI_Recv(y, input_array_size, MPI_INT, rank + m, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            x = Batcher_merge(x, y, current_array_size, input_array_size);
            current_array_size = current_array_size + input_array_size;
        }
        m = 2 * m;
    }
    if (rank == 0) {
       std::memcpy(container, x, array_size * sizeof(int));
    }

    return container;
}

std::vector<int> getRandomVector(int size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> vec(size);
    for (int i = 0; i < size; i++) { vec[i] = gen() % 100; }
    return vec;
}
