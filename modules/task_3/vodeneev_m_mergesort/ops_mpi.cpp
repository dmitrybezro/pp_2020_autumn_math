// Copyright 2020 Mikhail Vodeneev
#include <mpi.h>
#include <math.h>
#include <random>
#include <ctime>
#include <vector>
#include <algorithm>
#include "../../../modules/test_tasks/test_mpi/ops_mpi.h"

int compare(const void* x1, const void* x2) {
    return *(static_cast<const int*>(x1))
        - *(static_cast<const int*>(x2));
}

void quicksort(std::vector<int>* Arr, int first, int last) {
    int mid, tmp;
    int f = first, l = last;
    mid = (*Arr)[(f + l) / 2];
    do {
        while ((*Arr)[f] < mid)
            f++;
        while ((*Arr)[l] > mid)
            l--;
        if (f <= l) {
            tmp = (*Arr)[f];
            (*Arr)[f] = (*Arr)[l];
            (*Arr)[l] = tmp;
            f++;
            l--;
        }
    } while (f < l);
    if (first < l)
        quicksort(Arr, first, l);
    if (f < last)
        quicksort(Arr, f, last);
}

std::vector<int> merge(std::vector<int> a, std::vector<int> b) {
    int size1 = a.size();
    int size2 = b.size();
    std::vector<int> res(size1 + size2);
    int i = 0;
    int j = 0;
    int k = 0;
    while (i < size1 && j < size2) {
        if (a[i] <= b[j]) {
            res[k] = a[i];
            i++;
        } else {
            res[k] = b[j];
            j++;
        }
        k++;
    }
    while (i < size1) {
        res[k] = a[i];
        i++;
        k++;
    }

    while (j < size2) {
        res[k] = b[j];
        j++;
        k++;
    }
    return res;
}


std::vector<int> mergesort(std::vector<int> vec) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int delta = vec.size() / size;
    int delta_over = vec.size() % size;
    if (delta == 0) {
        std::vector<int> er(vec.size());
        size = vec.size();
        delta = 1;
        delta_over = 0;
        if (rank >= size)
            return er;
    }
    int delta_over2 = delta_over;
    if (rank != size - 1 || rank == 0)
        delta_over = 0;
    if (rank == 0) {
        for (int i = 1; i < size - 1; i++) {
            MPI_Send(&vec[0] + i * delta,
                    delta, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        if (size > 1)
        MPI_Send(&vec[0] + (size-1) * delta,
            delta + delta_over2, MPI_INT, size - 1, 0, MPI_COMM_WORLD);
    }
    std::vector<int> local_vec(delta + delta_over);
    if (rank == 0) {
        local_vec = std::vector<int>(vec.begin(), vec.begin() + delta);
        quicksort(&local_vec, 0, local_vec.size() - 1);
    } else {
        MPI_Status status;
        MPI_Recv(local_vec.data(), delta + delta_over,
            MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        quicksort(&local_vec, 0, local_vec.size() - 1);
    }

    int sizetmp = size;
    int k = 1;
    int n = delta + delta_over;

    while (sizetmp > 1) {
        sizetmp = sizetmp / 2 + sizetmp % 2;
        if ((rank - k) % (2 * k) == 0) {
            MPI_Send(&n, 1, MPI_INT, (rank - k), 0, MPI_COMM_WORLD);
            MPI_Send(local_vec.data(), n, MPI_INT, rank - k, 0, MPI_COMM_WORLD);
        }

        if ((rank % (2 * k) == 0) && (size - rank > k)) {
            MPI_Status status1, status2;
            int local_size;
            MPI_Recv(&local_size, 1, MPI_INT, rank + k,
                0, MPI_COMM_WORLD, &status1);
            std::vector <int> local_vec_2(local_size);
            MPI_Recv(local_vec_2.data(), local_size,
                MPI_INT, rank + k, 0, MPI_COMM_WORLD, &status2);
            local_vec = merge(local_vec, local_vec_2);
            n = n + local_size;
        }
        k = 2 * k;
    }
    return local_vec;
}
