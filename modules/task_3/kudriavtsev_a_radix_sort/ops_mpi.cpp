// Copyright 2020 Kudriavtsev Alexander
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <iostream>
#include "./ops_mpi.h"

std::vector<int> getRandomVector(int n, int mod) {  
    std::vector<int> vec(n);
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    for (int i = 0; i < n; ++i) {
        vec[i] = gen() % mod;  // Always is positive number. But tmp=gen(); tmp2=tmp%mod could be negative. MAGIC
    }
    return vec;
}

int* merge(int* vec1, int n1, int* vec2, int n2) {
    const int n = n1 + n2;
    int* vec = new int[n];
    int i = 0, j = 0, k = 0;
    while (j < n1 && k < n2) {
        if (vec1[j] < vec2[k])
            vec[i] = vec1[j++];
        else
            vec[i] = vec2[k++];
        ++i;
    }
    while (j < n1) {
        vec[i++] = vec1[j++];
    }
    while (k < n2) {
        vec[i++] = vec2[k++];
    }
    delete[] vec1;
    delete[] vec2;
    return vec;
}

int* parallelRadixSort(int* vec, int n)
{
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double t_b, t_e;
    if (rank == 0) {
        t_b = MPI_Wtime();
    }
    const int delta = n / size;
    const int rem = n % size;
    int* counts = new int[size * 2];  // Array of pairs would be faster(?) than this
    int* sendcounts = counts, * displs = counts + size;
    displs[0] = 0;
    sendcounts[0] = delta + static_cast<int>(0 < rem);
    for (int i = 1; i < size; ++i) {
        displs[i] = displs[i - 1] + sendcounts[i - 1];
        sendcounts[i] = delta + static_cast<int>(i < rem);
    }
    int proc_size = delta + static_cast<int>(rank < rem);
    int* arr_proc = new int[proc_size];
    MPI_Scatterv(vec, sendcounts, displs, MPI_INT, arr_proc, proc_size, MPI_INT, 0, MPI_COMM_WORLD);
    // -- Sorting --
    radixSort(arr_proc, proc_size);
    // -- Merging --
    int s = size, m = 1;
    while (s > 1) {
        s = s / 2 + s % 2;
        if ((rank - m) % (2 * m) == 0){
            MPI_Send(&proc_size, 1, MPI_INT,
                rank - m, 0, MPI_COMM_WORLD);
            MPI_Send(arr_proc, proc_size, MPI_INT,
                rank - m, 0, MPI_COMM_WORLD);
        }
        if ((rank % (2 * m) == 0) && ((size - rank) > m)){
            MPI_Status status;
            int k1;
            MPI_Recv(&k1, 1, MPI_INT,
                rank + m, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int* y = new int[k1];
            MPI_Recv(y, k1, MPI_INT,
                rank + m, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            arr_proc = merge(arr_proc, proc_size, y, k1);
            proc_size = proc_size + k1;
        }
        m = 2 * m;
    }
    if (rank == 0) {
        for (int i = 0; i < n; ++i) {
            vec[i] = arr_proc[i];
            //std::cout << arr_proc[i] << " ";
        }
        //std::cout << '\n';
        t_e = MPI_Wtime();
        std::cout << "Parallel time: " << t_e - t_b << std::endl;
    }
    delete[] arr_proc;
    return vec;
}

int* radixSort(int vec[], int n)
{
    int maxel = vec[0];
    const int power = 10;
    // --finding maximal element--
    for (int i = 0; i < n; ++i)
        if (vec[i] > maxel)
            maxel = vec[i];
    int* tmp_arr = new int[n];
    // -- Actual algoritm --
    for (int exp = 1; maxel / exp > 0; exp *= power) {
        int counts[power] = { 0 };  // if "power" is const, then static array would be more efficient
        for (int i = 0; i < n; ++i)
            ++counts[(vec[i] / exp) % power];  // Could use dynamic lists, but that (probably) would be s*cks
        for (int i = 1; i < power; ++i)
            counts[i] += counts[i-1];
        for (int i = n - 1; i >= 0; --i) {
            tmp_arr[counts[(vec[i] / exp) % power] - 1] = vec[i];
            --counts[(vec[i] / exp) % power];
        }
        for (int i = 0; i < n; ++i)
            vec[i] = tmp_arr[i];
    }
    // -- RELEASE RECOLLEC..HM! MEMORY --
    delete[] tmp_arr;
    return vec;
}

