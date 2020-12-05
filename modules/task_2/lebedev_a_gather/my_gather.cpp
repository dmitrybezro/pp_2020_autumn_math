// Copyright 2020 Lebedev Andrew
#include <math.h>
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <cstring>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_2/lebedev_a_gather/my_gather.h"

int my_gather(void* sbuf, int scount, MPI_Datatype stype, void* rbuf,
    int rcount, MPI_Datatype rtype, int dest, MPI_Comm comm) {
    int size, rank;
    if (rcount == 0) {
        return 0;
    }
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if (size == 1) {
        printf("Error: size must be > 1\n");
        return -1;
    }

    if (stype == MPI_INT) {
        std::memcpy(rbuf, sbuf, rcount * sizeof(int));
    } else if (stype == MPI_FLOAT) {
        std::memcpy(rbuf, sbuf, rcount * sizeof(float));
    } else if (stype == MPI_DOUBLE) {
        std::memcpy(rbuf, sbuf, rcount * sizeof(double));
    }
    MPI_Status status;
    int loc_count = rcount;
    int root_count = rcount;
    int curr_size = size;
    if (rtype == MPI_INT) {
        while (curr_size > 1) {
            int full = false;
            if (2 * rank + 2 < curr_size) {
                full = true;
            }
            if (rank == 0) {  // recieve only
                MPI_Recv(&loc_count, 1, MPI_INT, 1, 0, comm, &status);
                MPI_Recv(static_cast<int*>(rbuf) + root_count,
                    loc_count, rtype, 1, 0, comm, &status);
                root_count += loc_count;  // updating rbuf size on root

                if (curr_size != 2) {
                    MPI_Recv(&loc_count, 1, MPI_INT, 2, 0, comm, &status);
                    MPI_Recv(static_cast<int*>(rbuf) + root_count,
                        loc_count, rtype, 2, 0, comm, &status);
                    root_count += loc_count;  // updating rbuf size on root
                }

            } else if ((2 * (rank + 1)) > curr_size && rank < curr_size) {  // send only
                MPI_Send(&loc_count, 1, MPI_INT, (rank - 1) / 2 , 0, comm);
                MPI_Send(rbuf, loc_count, rtype, (rank - 1) / 2, 0, comm);

            } else if (rank < curr_size) {  // Send and Receive
                // Send Phase
                MPI_Send(&loc_count, 1, MPI_INT, (rank - 1) / 2 , 0, comm);
                MPI_Send(rbuf, loc_count, rtype, (rank - 1) / 2, 0, comm);

                // Receive Phase
                MPI_Recv(&loc_count, 1, MPI_INT, 2 * rank + 1, 0, comm, &status);
                MPI_Recv(static_cast<int*>(rbuf), loc_count, rtype,
                    2 * rank + 1, 0, comm, &status);
                if (full) {
                    int temp_count = loc_count;
                    MPI_Recv(&loc_count, 1, MPI_INT, 2 * rank + 2, 0, comm, &status);
                    MPI_Recv(static_cast<int*>(rbuf) + temp_count, loc_count, rtype,
                        2 * rank + 2, 0, comm, &status);
                    loc_count += temp_count;
                }
            }
            curr_size /= 2;
        }

        if (dest != 0) {
            if (rank == 0) {
                MPI_Send(rbuf, root_count, rtype, dest, 0, comm);
            } else if (rank == dest) {
                MPI_Recv(static_cast<int*>(rbuf),
                    size * rcount, rtype, 0, 0, comm, &status);
            }
        }
    } else if (rtype == MPI_FLOAT) {
        while (curr_size > 1) {
            int full = false;
            if (2 * rank + 2 < curr_size) {
                full = true;
            }
            /* Receive only */
            if (rank == 0) {
                MPI_Recv(&loc_count, 1, MPI_INT, 1, 0, comm, &status);
                MPI_Recv(static_cast<float*>(rbuf) + root_count,
                    loc_count, rtype, 1, 0, comm, &status);
                root_count += loc_count;  // updating rbuf size on root

                if (curr_size != 2) {
                    MPI_Recv(&loc_count, 1, MPI_INT, 2, 0, comm, &status);
                    MPI_Recv(static_cast<float*>(rbuf) + root_count,
                        loc_count, rtype, 2, 0, comm, &status);
                    root_count += loc_count;  // updating rbuf size on root
                }

            /* Send only */
            } else if ((2 * (rank + 1)) > curr_size && rank < curr_size) {
                MPI_Send(&loc_count, 1, MPI_INT, (rank - 1) / 2 , 0, comm);
                MPI_Send(rbuf, loc_count, rtype, (rank - 1) / 2, 0, comm);

            /* Send and receive */
            } else if (rank < curr_size) {
                // Send Phase
                MPI_Send(&loc_count, 1, MPI_INT, (rank - 1) / 2 , 0, comm);
                MPI_Send(rbuf, loc_count, rtype, (rank - 1) / 2, 0, comm);

                // Receive Phase
                MPI_Recv(&loc_count, 1, MPI_INT, 2 * rank + 1, 0, comm, &status);
                MPI_Recv(static_cast<float*>(rbuf), loc_count, rtype,
                    2 * rank + 1, 0, comm, &status);
                if (full) {
                    int temp_count = loc_count;
                    MPI_Recv(&loc_count, 1, MPI_INT, 2 * rank + 2, 0, comm, &status);
                    MPI_Recv(static_cast<float*>(rbuf) + temp_count, loc_count, rtype,
                        2 * rank + 2, 0, comm, &status);
                    loc_count += temp_count;
                }
            }
            curr_size /= 2;
        }

        if (dest != 0) {
            if (rank == 0) {
                MPI_Send(rbuf, root_count, rtype, dest, 0, comm);
            } else if (rank == dest) {
                MPI_Recv(static_cast<float*>(rbuf),
                    size * rcount, rtype, 0, 0, comm, &status);
            }
        }
    } else if (rtype == MPI_DOUBLE) {
                while (curr_size > 1) {
            int full = false;
            if (2 * rank + 2 < curr_size) {
                full = true;
            }
            if (rank == 0) {  // recieve only
                MPI_Recv(&loc_count, 1, MPI_INT, 1, 0, comm, &status);
                MPI_Recv(static_cast<double*>(rbuf) + root_count,
                    loc_count, rtype, 1, 0, comm, &status);
                root_count += loc_count;  // updating rbuf size on root

                if (curr_size != 2) {
                    MPI_Recv(&loc_count, 1, MPI_INT, 2, 0, comm, &status);
                    MPI_Recv(static_cast<double*>(rbuf) + root_count,
                        loc_count, rtype, 2, 0, comm, &status);
                    root_count += loc_count;  // updating rbuf size on root
                }

            } else if ((2 * (rank + 1)) > curr_size && rank < curr_size) {  // send only
                MPI_Send(&loc_count, 1, MPI_INT, (rank - 1) / 2 , 0, comm);
                MPI_Send(rbuf, loc_count, rtype, (rank - 1) / 2, 0, comm);

            } else if (rank < curr_size) {  // Send and Receive
                // Send Phase
                MPI_Send(&loc_count, 1, MPI_INT, (rank - 1) / 2 , 0, comm);
                MPI_Send(rbuf, loc_count, rtype, (rank - 1) / 2, 0, comm);

                // Receive Phase
                MPI_Recv(&loc_count, 1, MPI_INT, 2 * rank + 1, 0, comm, &status);
                MPI_Recv(static_cast<double*>(rbuf), loc_count, rtype,
                    2 * rank + 1, 0, comm, &status);
                if (full) {
                    int temp_count = loc_count;
                    MPI_Recv(&loc_count, 1, MPI_INT, 2 * rank + 2, 0, comm, &status);
                    MPI_Recv(static_cast<double*>(rbuf) + temp_count, loc_count, rtype,
                        2 * rank + 2, 0, comm, &status);
                    loc_count += temp_count;
                }
            }
            curr_size /= 2;
        }

        if (dest != 0) {
            if (rank == 0) {
                MPI_Send(rbuf, root_count, rtype, dest, 0, comm);
            } else if (rank == dest) {
                MPI_Recv(static_cast<double*>(rbuf),
                    size * rcount, rtype, 0, 0, comm, &status);
            }
        }
    } else {
        printf("Error: datatype not supported\n");
        return -1;
    }
    return MPI_SUCCESS;
}
