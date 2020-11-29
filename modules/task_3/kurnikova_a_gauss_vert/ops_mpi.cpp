// Copyright 2020 Kurnikova Anastasia
#include <mpi.h>
#include <vector>
#include "../../../modules/task_3/kurnikova_a_gauss_vert/ops_mpi.h"

std::vector<double> sequential(std::vector<double> image, int w, int h) {
    int it = 0;
    std::vector<double> res(w * h);
    std::vector<double> tmp((w + 2) * (h + 2));
    for (int x = 0; x < w + 2; x++)
        for (int y = 0; y < h + 2; y++)
            if ((x == 0) || (y == 0) || (x == w + 1) || (y == h + 1))
                tmp[x * (w + 2) + y] = 0;
            else
                tmp[x * (w + 2) + y] = image[(x - 1) * w + y - 1];
    for (int x = 1; x < w + 1; x++)
        for (int y = 1; y < w + 1; y++) {
            double sum = 0;
            for (int j = -1; j < 2; ++j)
                for (int k = -1; k < 2; ++k)
                    sum = sum + tmp[(x + j) * w + y + k] * gauss[j + k + 2];
            res[it] = sum / 16;
            ++it;
        }
    return res;
}

std::vector<double> parallel(std::vector<double> image, int n) {
    int procrank, procnum;
    MPI_Comm_size(MPI_COMM_WORLD, &procnum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procrank);
    int oneprocsize = n  / procnum;
    int oneprocdel = n % procnum;
    int firstprocdel = oneprocdel;
    if (procrank == 0) {
        for (int i = 1; i < procnum - 1; i++)
            MPI_Send(image.data() + i * oneprocsize, oneprocsize * n,
                                      MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        if (procnum > 1)
            MPI_Send(image.data() + (procnum - 1) * oneprocsize,
          (oneprocsize + oneprocdel) * n, MPI_DOUBLE, procnum - 1, 0,
                                                        MPI_COMM_WORLD);
    }
    if ((procrank != procnum - 1) || (procnum == 1))
        oneprocdel = 0;
    std::vector<double> oneprocim(n * (oneprocsize + oneprocdel));
    if (procrank == 0) {
        oneprocim = std::vector<double>(image.begin(),
                         image.begin() + oneprocsize * n);
        oneprocim = sequential(oneprocim, oneprocsize, n);
    } else {
        MPI_Status status;
        MPI_Recv(oneprocim.data(), n * (oneprocsize + oneprocdel),
                            MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        oneprocim = sequential(oneprocim, oneprocsize + oneprocdel, n);
    }
    std::vector<double> res(n * n);
    MPI_Gather(oneprocim.data(), oneprocsize * n, MPI_DOUBLE,
           res.data(), oneprocsize * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if ((procrank == procnum - 1) && (oneprocdel > 0))
        MPI_Send(oneprocim.data() + oneprocsize * n, oneprocdel * n,
                                     MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    if ((procrank == 0) && (oneprocdel > 0)) {
        MPI_Status status;
        MPI_Recv(res.data() + (procnum - 1) * oneprocsize * n,
                     firstprocdel * n, MPI_DOUBLE, procnum - 1, 0,
                          MPI_COMM_WORLD, &status);
        }
    return res;
}
