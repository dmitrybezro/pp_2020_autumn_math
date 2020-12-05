// Copyright 2020 Smirnov Aleksandr
#include <mpi.h>
#include <iostream>
#include <vector>
#include "../../../modules/task_2/smirnov_a_lattice_torus/lattice_torus.h"

MPI_Comm createTopology(MPI_Comm old_comm, int ndims) {
    int nnodes;
    MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
    std::vector<int> dims(ndims);
    MPI_Dims_create(nnodes, ndims, dims.data());
    std::vector<int> periods(ndims);
    periods.assign(ndims, 1);
    MPI_Comm cart_comm;
    MPI_Cart_create(old_comm, ndims, dims.data(), periods.data(), 0, &cart_comm);
    return cart_comm;
}

void Send(void* buf, int count, MPI_Datatype type, int* dest_coords, int tag, MPI_Comm comm) {
    int dest_rank, rank;
    MPI_Cart_rank(comm, dest_coords, &dest_rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == dest_rank)
        return;
    MPI_Send(buf, count, type, dest_rank, tag, comm);
}

void Recv(void* buf, int count, MPI_Datatype type, int* src_coords, int tag, MPI_Comm comm, MPI_Status* status) {
    int src_rank, rank;
    MPI_Cart_rank(comm, src_coords, &src_rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == src_rank)
        return;
    MPI_Recv(buf, count, type, src_rank, tag, comm, status);
}

void SendRecv(void* sbuf, int scount, MPI_Datatype stype, int* dest_coords, int stag,
    void* rbuf, int rcount, MPI_Datatype rtype, int* src_coords, int rtag, MPI_Comm comm, MPI_Status* status) {
    int dest_rank, src_rank, rank;
    MPI_Cart_rank(comm, dest_coords, &dest_rank);
    MPI_Cart_rank(comm, src_coords, &src_rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (dest_rank == src_rank)
        return;
    if (rank == src_rank) {
        MPI_Send(sbuf, scount, stype, dest_rank, stag, comm);
    }
    if (rank == dest_rank) {
        MPI_Recv(rbuf, rcount, rtype, src_rank, rtag, comm, status);
    }
}
