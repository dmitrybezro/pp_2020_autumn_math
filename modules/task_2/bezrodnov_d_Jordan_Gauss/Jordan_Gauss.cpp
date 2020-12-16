// Copyright 2020 Bezrodnov Dmitry
#include <math.h>
#include<mpi.h>
#include<iostream>
#include<vector>
#include <ctime>
#include <cmath>
#include <limits>
#include <utility>
#include <random>
#include "../../../modules/task_2/bezrodnov_d_Jordan_Gauss/Jordan_Gauss.h"

std::vector<double> getRandomMatrix(int size) {
    // диапазон от min до max
    double min = -1488.415;
    double max = 988.947;
    std::vector<double> matr(size * (size + 1));
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));

    for (size_t i = 0; i < matr.size(); i++) {
        matr[i] = static_cast<double>(gen()) / (1000 * (max - min)) + min * pow(-1, i);
    }

    return matr;
}

std::vector<double> ParallelJordanGauss(const std::vector<double>& MainMatr, int size) {
    int RANK;
    int SIZE;
    MPI_Comm_size(MPI_COMM_WORLD, &SIZE);
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);

    if (SIZE == 1) {
        return SequenJordanGauss(MainMatr, size);
    }

    const int div = size / (SIZE - 1);
    const int mod = size % (SIZE - 1);

    // Столбцовые части большого вектора
    std::vector<double> part_vec(div * size);

    std::vector<double> koef_recv;
    std::vector<double> koef_send;
    std::vector<double> main_elem(size);
    std::vector<double> X(size);
    double answer;

    if (RANK == 0) {
        for (int proc = 0; proc < SIZE - 2; proc++) {
            MPI_Send(&MainMatr[0] + div * size * proc, div * size, MPI_DOUBLE, proc + 1, 0, MPI_COMM_WORLD);
        }
        // 1 + div потому что еще столбец b
        MPI_Send(&MainMatr[0] + MainMatr.size() - size * (mod + div + 1),
            size * (mod + div + 1), MPI_DOUBLE, SIZE - 1, 0, MPI_COMM_WORLD);
    } else {
        MPI_Status stat;
        if (RANK < SIZE - 1) {
            MPI_Recv(&part_vec[0], div * size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);
        } else {
            part_vec.resize(size * (div + mod + 1));
            MPI_Recv(&part_vec[0], (div + mod + 1) * size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);
        }
    }

    //  Начинаем вычисления
    if (RANK != 0) {
        //  Количество исходных столбцов в векторе
        int number_columns;
        number_columns = part_vec.size() / size;

        koef_recv.resize(size);
        koef_send.resize(size);
        int pos_main_elem = 0;

        //  Принять приходящие коэфы и с ними обратботать вектора
        //  В первый процесс не откуда принимать коэфы
            for (int proc = 0; proc < RANK - 1; proc++) {
                //  Приняли коэфы от других процессов
                for (int num_col_koef = 0; num_col_koef < div; num_col_koef++) {
                    MPI_Status stat;
                    MPI_Recv(&koef_recv[0], size, MPI_DOUBLE, proc + 1, 0, MPI_COMM_WORLD, &stat);

                    //  Нашли позицию главного элемента в пришедшем векторе коэффициентов
                    for (size_t i = 0; i < koef_recv.size(); i++) {
                        if (koef_recv[i] == 999) {
                            pos_main_elem = i;
                        }
                    }

                    //  Обработали вектор процесса с пришедшими коэфами
                    for (int j = 0; j < number_columns; j++) {
                        for (int i = 0; i < size; i++) {
                            if (i != pos_main_elem) {
                                part_vec[i + j * size] -= koef_recv[i] * part_vec[pos_main_elem + j * size];
                            }
                        }
                    }
                }
            }

        int vec;
        for (vec = 0; vec < number_columns; vec++) {
            if (((RANK == (SIZE - 1)) && (vec == (div + mod)))) {
                break;
            }
            if ((RANK != 1) || (vec != 0)) {
                pos_main_elem++;
            }

            //  Заполнили один столбец вектора коэффициентами
            for (int i = 0; i < size; i++) {
                if (i == pos_main_elem) {
                    koef_send[i] = 999;  //  ДЛя всех кроме первого
                } else {
                    koef_send[i] = part_vec[i + vec * size] / part_vec[pos_main_elem + vec * size];
                }
            }

            //  Привели все столбцы в векторе процесса, к нужному виду( ноль под ведущим элементом)
            for (int j = vec; j < number_columns; j++) {
                for (int i = 0; i < size; i++) {
                    if (i != pos_main_elem) {
                        if ((i > pos_main_elem) && (j == vec)) {
                            part_vec[i + j * size] = 0;
                        } else {
                            part_vec[i + j * size] -= koef_send[i] * part_vec[pos_main_elem + j * size];
                        }
                    }
                }
            }

            //  В последнем процессе отправки делать не нужно
            if (RANK != (SIZE - 1)) {
                //  Отправили каждый главный элемент вектора в последний процесс
                answer = part_vec[pos_main_elem + vec * size];
                MPI_Send(&answer, 1, MPI_DOUBLE, SIZE - 1, 1, MPI_COMM_WORLD);

                //  Раздать коэффициенты столбцов остальным процессам
                for (int proc = RANK + 1; proc < SIZE; proc++) {
                    MPI_Send(&koef_send[0], koef_send.size(), MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
                }
            } else {
                main_elem[main_elem.size() - (number_columns - 1) + vec] = part_vec[pos_main_elem + vec * size];
            }
        }
        //  Обработка ответа
        if ((RANK == (SIZE - 1)) && (vec == (div + mod))) {
            int m = 0;
            for (int i = 0; i < SIZE - 2; i++) {
                //  Приняли остальные диагональные элементы
                for (int num_col_koef = 0; num_col_koef < div; num_col_koef++, m++) {
                    MPI_Status stat;
                    MPI_Recv(&main_elem[m], 1, MPI_DOUBLE, i + 1, 1, MPI_COMM_WORLD, &stat);
                }
            }
            //  Поделили полученный столбец b, на диагональные элементы
            for (int j = size, i = 0; j > 0; j--, i++) {
                part_vec[part_vec.size() - j] /= main_elem[i];
                X[i] = part_vec[part_vec.size() - j];
            }
            MPI_Send(&X[0], X.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
    if (RANK == 0) {
        MPI_Status stat;
        MPI_Recv(&X[0], size, MPI_DOUBLE, SIZE - 1, 0, MPI_COMM_WORLD, &stat);
    }
    return X;
}

std::vector<double> SequenJordanGauss(const std::vector<double>& MainMatr, int size) {
    std::vector<double> koef_send(size);
    std::vector<double> X(size);
    std::vector<double> part_vec(MainMatr);
    for (int vec = 0; vec < size; vec++) {
        int pos_main_elem = vec;
        //  Заполнили один столбец вектора коэффициентами
        for (int i = 0 ; i < size; i++) {
            if (i == pos_main_elem) {
                koef_send[i] = 0;  //  ДЛя всех кроме первого
            } else {
                koef_send[i] = part_vec[i + vec * size] / part_vec[pos_main_elem + vec * size];
            }
        }
        //  Привели все столбцы в векторе процесса, к нужному виду( ноль под ведущим элементом)
        for (int j = 0; j < size + 1; j++) {
            for (int i = 0; i < size; i++) {
                if (i != pos_main_elem) {
                    if ((i > pos_main_elem) && (j == vec)) {
                        part_vec[i + j * size] = 0;
                    } else {
                        part_vec[i + j * size] -= koef_send[i] * part_vec[pos_main_elem + j * size];
                    }
                }
            }
        }
    }

    for (int i = size, j = 0, m = 0; i > 0; i--, j +=(size + 1), m++) {
        part_vec[part_vec.size() - i]/=part_vec[j];
        X[m] = part_vec[part_vec.size() - i];
    }

    return X;
}

std::vector<double> MultiInverseMatrix(const std::vector<double>& MainMatr, int size) {
    std::vector<double> MatrWithout(MainMatr.size() - size);
    std::vector<double> b(size);

    for (size_t i = 0; i < MatrWithout.size(); i++) {
        MatrWithout[i] = MainMatr[i];
    }

    for (int i = size, j = 0; i > 0; i--, j++) {
        b[j] = MainMatr[MainMatr.size() - i];
    }

    double **matr = new double*[size];
    matr = VectorInPointer(MatrWithout, size);

    double** matrT = new double*[size];
    for (int i = 0; i < size; i++) {
        matrT[i] = new double[size];
    }

    Transpon(matr, matrT, size);

    std::vector<double> InvMatr = Inverse(matrT, size);
    std::vector<double> X = MultiMatrVector(InvMatr, b);
    return X;
}

void FreeMem(double **matr, int n) {
    for (int i = 0; i < n; i++)
        delete [] matr[i];
    delete [] matr;
}

//  функция вычеркивания строки и столбца
void Get_matr(double **matr, int n, double **temp_matr, int indRow, int indCol) {
    int ki = 0;
    for (int i = 0; i < n; i++) {
        if (i != indRow) {
            for (int j = 0, kj = 0; j < n; j++) {
                if (j != indCol) {
                    temp_matr[ki][kj] = matr[i][j];
                    kj++;
                }
            }
            ki++;
       }
    }
}

//  функция вычисления определителя матрицы
double Det(double **matr, int n) {
    double temp = 0;   //  временная переменная для хранения определителя
    int k = 1;      //  степень
    if (n == 1) {
        temp = matr[0][0];
    } else {
        if (n == 2) {
            temp = matr[0][0] * matr[1][1] - matr[1][0] * matr[0][1];
        } else {
            for (int i = 0; i < n; i++) {
                int m = n - 1;
                double **temp_matr = new double*[m];
                for (int j = 0; j < m; j++)
                    temp_matr[j] = new double[m];
                Get_matr(matr, n, temp_matr, 0, i);
                temp = temp + k * matr[0][i] * Det(temp_matr, m);
                k = -k;
                FreeMem(temp_matr, m);
            }
        }
    }
    return temp;
}

//  функция переводит длинный вектор в указатель
double** VectorInPointer(const std::vector<double> _matr, int cols) {
        double **matr = new double * [cols];
        for (int i = 0; i < cols; i++) {
            matr[i] = new double[cols];
        }

        for (int i = 0; i < _matr.size()/cols; i++) {
            for (int j = 0; j < cols; j++) {
                matr[i][j] = _matr[j + i*cols];
            }
        }
        return matr;
}

void Transpon(double **matr, double **tMatr, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            tMatr[j][i] = matr[i][j];
}


std::vector<double> Inverse(double** matr, int n) {
    double **obr_matr = new double * [n];
    double **tobr_matr = new double * [n];
    for (int i = 0; i < n; i++) {
        obr_matr[i] = new double[n];
        tobr_matr[i] = new double[n];
    }

    double det = Det(matr, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int m = n - 1;
                double **temp_matr = new double * [m];
                for (int k = 0; k < m; k++)
                    temp_matr[k] = new double[m];
                Get_matr(matr, n, temp_matr, i, j);
                obr_matr[i][j] = pow(-1.0, i + j + 2) * Det(temp_matr, m) / det;
                FreeMem(temp_matr, m);
            }
        }
    Transpon(obr_matr, tobr_matr, n);
    std::vector<double> res;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res.push_back(tobr_matr[i][j]);
        }
    }
    FreeMem(obr_matr, n);
    FreeMem(tobr_matr, n);
    return res;
}

std::vector<double> MultiMatrVector(const std::vector<double>& matr, const std::vector<double>& vec) {
    std::vector<double> result;
    result.assign(vec.size(), 0);
    for (size_t i = 0; i < vec.size(); i++) {
        for (int j = 0; j < vec.size(); j++) {
            result[i] += matr[j + vec.size() * i] * vec[j];
        }
    }
    return result;
}
