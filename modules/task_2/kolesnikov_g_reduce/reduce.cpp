// Copyright 2020 Kolesnikov Gleb
#include <mpi.h>
#include <string>
#include <ctime>
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include "../../../modules/task_2/kolesnikov_g_reduce/reduce.h"
int MPI_Reduce_My_Own(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
    MPI_Op op, int root, MPI_Comm comm) {
    int size;
    int rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    if (datatype == MPI_INT) {
        memcpy(recvbuf, sendbuf, count * sizeof(int));
    }
    int level = log2f(size + 1);
    int cur_n = size;
    int delta = exp2(level - 1);
    if (cur_n != 2 * delta - 1) {
        if (rank >= 2 * delta - 1 && rank < cur_n) {
            MPI_Send(recvbuf, count, datatype, 2 * (2 * delta - 1) - rank - 1, 0, comm);
        }
        if (rank >= 4 * delta - 2 - cur_n && rank < 2 * delta - 1) {
            MPI_Status status;
            if (datatype == MPI_INT) {
                int * temp_buf = new int[count];
                MPI_Recv(temp_buf, count, datatype, 4 * delta - 3 - rank, 0, comm, &status);
                if (op == MPI_SUM) {
                    for (int i = 0; i < count; i++) {
                        static_cast<int*>(recvbuf)[i] += temp_buf[i];
                    }
                } else {
                    if (op == MPI_MIN) {
                        for (int i = 0; i < count; i++) {
                            if (temp_buf[i] < static_cast<int*>(recvbuf)[i]) {
                                static_cast<int*>(recvbuf)[i] = temp_buf[i];
                            }
                        }
                    } else {
                        if (op == MPI_PROD) {
                            for (int i = 0; i < count; i++) {
                                static_cast<int*>(recvbuf)[i] *= temp_buf[i];
                            }
                        } else {
                            if (op == MPI_MAX) {
                                for (int i = 0; i < count; i++) {
                                    if (temp_buf[i] > static_cast<int*>(recvbuf)[i]) {
                                        static_cast<int*>(recvbuf)[i] = temp_buf[i];
                                    }
                                }
                            } else {
                                return -1;
                            }
                        }
                    }
                }
            }
        }
    }

    cur_n = 2 * delta - 1;

    while (level != 1) {
        if (rank >= cur_n / 2 && rank < cur_n) {
            MPI_Send(recvbuf, count, datatype, (rank - 1) / 2, 0, comm);
        }
        cur_n = cur_n - delta;
        delta = delta / 2;
        level--;
        if (rank >= cur_n - delta && rank < cur_n) {
            MPI_Status status;
            if (datatype == MPI_INT) {
                int * temp_buf = new int[count];
                if (op == MPI_SUM) {
                        MPI_Recv(temp_buf, count, datatype, 2 * rank + 1, 0, comm, &status);
                        for (int i = 0; i < count; i++) {
                            static_cast<int*>(recvbuf)[i] += temp_buf[i];
                        }
                        MPI_Recv(temp_buf, count, datatype, 2 * rank + 2, 0, comm, &status);
                        for (int i = 0; i < count; i++) {
                            static_cast<int*>(recvbuf)[i] += temp_buf[i];
                        }
                } else {
                    if (op == MPI_MIN) {
                        MPI_Recv(temp_buf, count, datatype, 2 * rank + 1, 0, comm, &status);
                        for (int i = 0; i < count; i++) {
                            if (temp_buf[i] < static_cast<int*>(recvbuf)[i]) {
                                static_cast<int*>(recvbuf)[i] = temp_buf[i];
                            }
                        }
                        MPI_Recv(temp_buf, count, datatype, 2 * rank + 2, 0, comm, &status);
                        for (int i = 0; i < count; i++) {
                            if (temp_buf[i] < static_cast<int*>(recvbuf)[i]) {
                                static_cast<int*>(recvbuf)[i] = temp_buf[i];
                            }
                        }
                    } else {
                        if (op == MPI_MAX) {
                            MPI_Recv(temp_buf, count, datatype, 2 * rank + 1, 0, comm, &status);
                            for (int i = 0; i < count; i++) {
                                if (temp_buf[i] > static_cast<int*>(recvbuf)[i]) {
                                    static_cast<int*>(recvbuf)[i] = temp_buf[i];
                                }
                            }
                            MPI_Recv(temp_buf, count, datatype, 2 * rank + 2, 0, comm, &status);
                            for (int i = 0; i < count; i++) {
                                if (temp_buf[i] > static_cast<int*>(recvbuf)[i]) {
                                    static_cast<int*>(recvbuf)[i] = temp_buf[i];
                                }
                            }
                        } else {
                            if (op == MPI_PROD) {
                                    MPI_Recv(temp_buf, count, datatype, 2 * rank + 1, 0, comm, &status);
                                    for (int i = 0; i < count; i++) {
                                        static_cast<int*>(recvbuf)[i] *= temp_buf[i];
                                    }
                                    MPI_Recv(temp_buf, count, datatype, 2 * rank + 2, 0, comm, &status);
                                    for (int i = 0; i < count; i++) {
                                        static_cast<int*>(recvbuf)[i] *= temp_buf[i];
                                    }
                                }
                        }
                    }
                }
            }
        }
    }

    if (rank == 0 && root != 0) {
        MPI_Send(recvbuf, count, datatype, root, 0, comm);
    }
    if (root != 0 && rank == root) {
        MPI_Status status;
        MPI_Recv(recvbuf, count, datatype, 0, 0, comm, &status);
    }

    return MPI_SUCCESS;
}
