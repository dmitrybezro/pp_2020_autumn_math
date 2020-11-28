// Copyright 2020 Smirnov Alexander
#ifndef  MODULES_TASK_2_SMIRNOV_A_LATTICE_TORUS_LATTICE_TORUS_H_
#define  MODULES_TASK_2_SMIRNOV_A_LATTICE_TORUS_LATTICE_TORUS_H_

MPI_Comm createTopology(MPI_Comm old_comm, int ndims);
void Send(void* buf, int count, MPI_Datatype type, int* dest_coords, int tag, MPI_Comm comm);
void Recv(void* buf, int count, MPI_Datatype type, int* src_coords, int tag, MPI_Comm comm, MPI_Status* status);
void SendRecv(void* sbuf, int scount, MPI_Datatype stype, int* dest_coords, int stag,
              void* rbuf, int rcount, MPI_Datatype rtype, int* src_coords, int rtag, MPI_Comm comm, MPI_Status* status);

#endif  //  MODULES_TASK_2_SMIRNOV_A_LATTICE_TORUS_LATTICE_TORUS_H_
