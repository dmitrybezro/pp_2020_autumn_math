// Copyright 2018 Smirnov Alexander
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <utility>
#include <tuple>
#include "../../../modules/task_3/smirnov_a_sparse_matrix_multiplication_ccs/sparse_matrix_multiplication_ccs.h"


std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> getRandomMatrix(int sz) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));

    int count = gen() % (sz * sz / 2 - 1);
    std::vector<double> values(count);
    std::vector<int> cols(count);
    std::vector<int> pointer(sz + 1, 0);
    for (int  i = 0; i < count; i++) { values[i] = gen() % 100; }
    for (int i = 0; i < count; i++) { cols[i] = gen() % sz; }
    for (int i = 1; i < sz + 1; i++) {
        if (i == sz) pointer[i] = count;
        while (true) {
            pointer[i] = pointer[i - 1] + gen() % (count / 2);
            if (pointer[i] <= count) break;
        }
    }
    return {values, cols, pointer};
}

std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> transpose(std::vector<double> values,
                                                                              std::vector<int> cols,
                                                                              std::vector<int> pointers) {
    int size = pointers.size() - 1;
    int count = values.size();

    std::vector<int> row;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < pointers[i + 1] - pointers[i]; j++)
            row.push_back(i);

    std::vector<double> valuesT;
    std::vector<int> colsT;
    std::vector<int> pointersT;

    pointersT.push_back(0);
    for (int i = 0; i < size; i++) {
        int point = 0;
        for (int j = 0; j < count; j++) {
            if (cols[j] == i) {
                valuesT.push_back(values[j]);
                colsT.push_back(row[j]);
                point++;
            }
        }
        pointersT.push_back(pointersT[pointersT.size() - 1] + point);
    }
    return { valuesT, colsT, pointersT };
}

std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
getSequentialOperations(std::vector<double> Avalues, std::vector<int> Acols, std::vector<int> Apointers,
    std::vector<double> Bvalues, std::vector<int> Bcols, std::vector<int> Bpointers) {
    int size = Apointers.size() - 1;
    std::vector<double> values;
    std::vector<int> cols;
    std::vector<int> pointers = { 0 };
    for (int i = 0; i < size; i++) {
        int coun = 0;
        for (int k = 0; k < size; k++) {
            std::vector<std::pair<double, int>> A, B;
            for (int j = 0; j < Apointers[i + 1] - Apointers[i]; j++) A.push_back(std::pair<double, int>
                                                                                  (Avalues[j + Apointers[i]],
                                                                                   Acols[j + Apointers[i]]));
            for (int j = 0; j < Bpointers[k + 1] - Bpointers[k]; j++) B.push_back(std::pair<double, int>
                                                                                  (Bvalues[j + Bpointers[k]],
                                                                                   Bcols[j + Bpointers[k]]));
            if (A.size() > 0 && B.size() > 0) {
                double sum = 0;
                int flag = 0;
                for (int j = 0; j < A.size(); j++)
                    for (int l = 0; l < B.size(); l++) {
                        if (A[j].second == B[l].second) {
                            flag++;
                            sum += A[j].first * B[l].first;
                        }
                    }
                if (flag) {
                    values.push_back(sum);
                    cols.push_back(k);
                    coun++;
                }
            }
        }
        if (coun) {
            pointers.push_back(pointers[pointers.size() - 1] + coun);
        } else {
            pointers.push_back(pointers[pointers.size() - 1]);
        }
    }
    return { values, cols, pointers };
}

std::vector<double> helpSeq4Par(std::vector<double> Avalues, std::vector<int> Acols,
                                std::vector<double> Bvalues, std::vector<int> Bcols) {
    std::vector<double> values;
    if (Avalues.size() > 0  && Bvalues.size() > 0) {
        double sum = 0.0;
        int flag = 0;
        for (int j = 0; j < Avalues.size(); j++)
            for (int l = 0; l < Bvalues.size(); l++) {
                if (Acols[j] == Bcols[l]) {
                    flag++;
                    sum += Avalues[j] * Bvalues[l];
                }
            }
        if (flag)
            values.push_back(sum);
    }
    return values;
}

std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
getParallelOperations(int sz, std::vector<double> Avalues, std::vector<int> Acols, std::vector<int> Apointers,
                              std::vector<double> Bvalues, std::vector<int> Bcols, std::vector<int> Bpointers) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<double> values;
    std::vector<int> cols;
    std::vector<int> pointers = { 0 };

    if (size < sz)
        return getSequentialOperations(Avalues, Acols, Apointers, Bvalues, Bcols, Bpointers);

    for (int i = 0; i < sz; i++) {
        int count = 0;
        std::vector<double> tmp;
        for (int k = 0; k < sz; k++) {
            if (k == rank) {
                tmp = helpSeq4Par(std::vector<double>(Avalues.begin()+Apointers[i], Avalues.begin()+Apointers[i + 1]),
                                  std::vector<int>(Acols.begin() + Apointers[i], Acols.begin() + Apointers[i + 1]),
                                  std::vector<double>(Bvalues.begin()+Bpointers[k], Bvalues.begin()+Bpointers[k + 1]),
                                  std::vector<int>(Bcols.begin() + Bpointers[k], Bcols.begin() + Bpointers[k + 1]));
            }
        }
        // std::cout << "WOW" << std::endl;
        if (rank == 0) {
            values.insert(values.end(), tmp.begin(), tmp.end());
            if (tmp.size() > 0) {
                cols.push_back(0);
                count += tmp.size();
            }
            for (int k = 1; k < sz; k++) {
                int tmp_size = 0;
                MPI_Status stat;
                // std::cout << "WOW - 0" << rank << std::endl;
                MPI_Recv(&tmp_size, 1, MPI_INT, k, 0, MPI_COMM_WORLD, &stat);
                if (tmp_size > 0) {
                    std::vector<double> tmp(tmp_size);
                    MPI_Recv(tmp.data(), tmp_size, MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &stat);
                    values.insert(values.end(), tmp.begin(), tmp.end());
                    cols.push_back(k);
                    count += tmp.size();
                }
            }
            if (count) {
                pointers.push_back(pointers[pointers.size() - 1] + count);
            } else {
                pointers.push_back(pointers[pointers.size() - 1]);
            }
        } else {
            if (rank < sz) {
                int tmp_size = tmp.size();
                MPI_Send(&tmp_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                if (tmp_size > 0)
                    MPI_Send(tmp.data(), tmp_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            }
        }
    }
    return { values, cols, pointers };
}
