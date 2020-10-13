// Copyright 2018 Nesterov Alexander
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include "../../../modules/test_tasks/test_mpi/ops_mpi.h"


std::vector<int> getRandomVector(int sz) {
	std::mt19937 gen;
	gen.seed(static_cast<unsigned int>(time(0)));
	std::vector<int> vec(sz);
	for (int i = 0; i < sz; i++) { vec[i] = gen() % 100 - 50; }
	return vec;
}

int getSequentialOperations(std::vector<int> vec) {
	const int  sz = vec.size();
	int change_of_sings = 0;
	bool positive;
	int i = 0;

	while (vec[i] == 0) i++;

	if (vec[i] > 0) positive = true;
	else positive = false;
	
	for (; i < sz; i++) {
		if (positive == true && vec[i] < 0) {
			positive = false;
			change_of_sings++;
		}
		if (positive == false && vec[i] > 0) {
			positive = true;
			change_of_sings++;
		}
	}
	return change_of_sings;
}

int getParallelOperations(std::vector<int> global_vec, int count_size_vector) {
	int size, rank, done = 0, * displs, * scounts;;
	double t1, t2;
	int n;
	int change_of_sings = 0, mychange = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Bcast(&global_vec, count_size_vector, MPI_INT, 0, MPI_COMM_WORLD);
	int local_count = 0;
	for (int i = rank; i < (count_size_vector - 1); i += size)
		local_count += getSequentialOperations(std::vector<int>{global_vec[i], global_vec[i + 1]});
	MPI_Reduce(&local_count, &change_of_sings, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	return change_of_sings;
}
