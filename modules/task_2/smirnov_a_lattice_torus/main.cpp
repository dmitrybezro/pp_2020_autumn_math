// Copyright 2020 Smirnov Alexander
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <ctime>
#include "./lattice_torus.h"

TEST(Parallel_Operations_MPI, Test_create_topology) {
    MPI_Comm comm_cart = createTopology(MPI_COMM_WORLD, 2);
    int status, ndims;
    MPI_Topo_test(comm_cart, &status);
    ASSERT_EQ(MPI_CART, status);

    MPI_Cartdim_get(comm_cart, &ndims);
    ASSERT_EQ(2, ndims);

    std::vector<int> dims(ndims), periods(ndims), coords(ndims);
    MPI_Cart_get(comm_cart, ndims, dims.data(), periods.data(), coords.data());
    EXPECT_TRUE(periods[0] & periods[1]);
}

TEST(Parallel_Operations_MPI, Test_send_from1_toLast) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm comm_cart = createTopology(MPI_COMM_WORLD, 2);

    std::vector<int> src_coords(2), dest_coords(2);
    MPI_Cart_coords(comm_cart, 0, 2, src_coords.data());
    MPI_Cart_coords(comm_cart, size - 1, 2, dest_coords.data());

    std::vector<int> buf(4);
    if (rank == 0) {
        buf = { 35, 876, 2, 78 };
        Send(buf.data(), 4, MPI_INT, dest_coords.data(), 0, comm_cart);
    }
    if (rank == size - 1) {
        MPI_Status status;
        Recv(buf.data(), 4, MPI_INT, src_coords.data(), 0, comm_cart, &status);
        ASSERT_EQ((std::vector<int>{ 35, 876, 2, 78 }), buf);
    }
}

TEST(Parallel_Operations_MPI, Test_sendrecv_from1_toLast) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm comm_cart = createTopology(MPI_COMM_WORLD, 2);

    std::vector<int> src_coords(2), dest_coords(2);
    MPI_Cart_coords(comm_cart, 0, 2, src_coords.data());
    MPI_Cart_coords(comm_cart, size - 1, 2, dest_coords.data());

    std::vector<int> buf(10);
    if (rank == 0) {
        buf = { 11, 22, 33, 44, 55, 66, 77, 88, 99, 100 };
    }
    MPI_Status status;
    SendRecv(buf.data(), 10, MPI_INT, dest_coords.data(), 0,
             buf.data(), 10, MPI_INT, src_coords.data(), 0, comm_cart, &status);
    if (rank == size - 1) {
        ASSERT_EQ((std::vector<int>{ 11, 22, 33, 44, 55, 66, 77, 88, 99, 100 }), buf);
    }
}

TEST(Parallel_Operations_MPI, Test_send_fromA_toB) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm comm_cart = createTopology(MPI_COMM_WORLD, 2);

    std::mt19937 generator;
    generator.seed(static_cast<unsigned int>(time(0)));
    int src_rank = generator() % size;
    int dest_rank = generator() % size;

    std::vector<int> src_coords(2), dest_coords(2);
    MPI_Cart_coords(comm_cart, src_rank, 2, src_coords.data());
    MPI_Cart_coords(comm_cart, dest_rank, 2, dest_coords.data());

    std::vector<int> buf(4);
    if (rank == src_rank) {
        buf = { 35, 876, 2, 78 };
        Send(buf.data(), 4, MPI_INT, dest_coords.data(), 0, comm_cart);
    }
    if (rank == dest_rank) {
        MPI_Status status;
        Recv(buf.data(), 4, MPI_INT, src_coords.data(), 0, comm_cart, &status);
        ASSERT_EQ((std::vector<int>{ 35, 876, 2, 78 }), buf);
    }
}

TEST(Parallel_Operations_MPI, Test_sendrecv_fromA_toB) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm comm_cart = createTopology(MPI_COMM_WORLD, 2);

    std::mt19937 generator;
    generator.seed(static_cast<unsigned int>(time(0)));
    int src_rank = generator() % size;
    int dest_rank = generator() % size;

    std::vector<int> src_coords(2), dest_coords(2);
    MPI_Cart_coords(comm_cart, src_rank, 2, src_coords.data());
    MPI_Cart_coords(comm_cart, dest_rank, 2, dest_coords.data());

    std::vector<int> buf(10);
    if (rank == src_rank) {
        buf = { 11, 22, 33, 44, 55, 66, 77, 88, 99, 100 };
    }
    MPI_Status status;
    SendRecv(buf.data(), 10, MPI_INT, dest_coords.data(), 0,
             buf.data(), 10, MPI_INT, src_coords.data(), 0, comm_cart, &status);
    if (rank == dest_rank) {
        ASSERT_EQ((std::vector<int>{ 11, 22, 33, 44, 55, 66, 77, 88, 99, 100 }), buf);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
