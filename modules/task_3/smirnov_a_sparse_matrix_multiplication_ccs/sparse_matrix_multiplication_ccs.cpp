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

std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
getParallelOperations(int sz, std::vector<double> Avalues, std::vector<int> Acols, std::vector<int> Apointers,
    std::vector<double> Bvalues, std::vector<int> Bcols, std::vector<int> Bpointers) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<double> values;
    std::vector<int> cols;
    std::vector<int> pointers = { 0 };
    std::vector<int> size_v(size);
    std::vector<int> size_c(size);
    std::vector<int> size_p(size);

    if (size < sz)
        return getSequentialOperations(Avalues, Acols, Apointers, Bvalues, Bcols, Bpointers);

    for (int i = 0; i < sz; i++) {
        if (i == rank) {
            int coun = 0;
            std::vector<double> local_values;
            std::vector<int> local_cols;
            std::vector<int> local_pointers;
            for (int k = 0; k < sz; k++) {
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
                        local_values.push_back(sum);
                        local_cols.push_back(k);
                        coun++;
                    }
                }
            }
            if (coun) {
                local_pointers.push_back(local_pointers[local_pointers.size()] + coun);
            } else {
                local_pointers.push_back(local_pointers[local_pointers.size()]);
            }
            int slv = local_values.size();
            int slc = local_cols.size();
            int slp = local_pointers.size();
            MPI_Gather(&slv, 1, MPI_INT, size_v.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(&slc, 1, MPI_INT, size_c.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(&slp, 1, MPI_INT, size_p.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
            std::vector<int> displs_c(size, 0);
            std::vector<int> displs_v(size, 0);
            std::vector<int> displs_p(size, 0);
            for (int j = 0, v = 0, c = 0, p = 0; j < size; j++) {
                v += size_v[j];
                c += size_c[j];
                p += size_p[j];
                displs_c[j] = c;
                displs_v[j] = v;
                displs_p[j] = p;
            }

            MPI_Gatherv(local_values.data(), local_values.size(), MPI_DOUBLE, values.data(),
                size_v.data(), displs_v.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gatherv(local_cols.data(), local_cols.size(), MPI_INT, cols.data(),
                size_c.data(), displs_c.data(), MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gatherv(local_pointers.data(), local_pointers.size(), MPI_INT, pointers.data(),
                size_p.data(), displs_p.data(), MPI_INT, 0, MPI_COMM_WORLD);
        }
    }
    return { values, cols, pointers };
}
