#include<mpi.h>
#include <iostream>
#include<vector>
#include <math.h>
#include <random>
#include <ctime>
#include <stdlib.h>

using namespace std;

struct point {
    double x;
    double y;
};

//  С какой стороны находит точка B относительно вектора A1A2
//  Или направление поворота: > 0 - поворот левый 
//  < 0 - поворот правый, = 0 - на одной прямой (коллинеарны)
double SideSpace(point A1, point A2, point B)
{
    return (A2.x - A1.x) * (B.y - A2.y) - (A2.y - A1.y) * (B.x - A2.x);
}

void InsertionSort(const vector<point>& _cloud, vector<int>& list_point)
{
    for (int i = 0; i < _cloud.size(); i++) {
        int j = i;
        while (j > 1 && SideSpace(_cloud[list_point[0]], _cloud[list_point[j - 1]], _cloud[list_point[j]]) < 0) {
            swap(list_point[j], list_point[j - 1]);
            j--;
        }
    }
}

vector<point> getRandomCloud(int size) {
    // диапазон от min до max
    double min = -1488.415;
    double max = 988.947;
    vector<point> cloud(size);
    srand(time(NULL));
    for (int i = 0; i < size; i++) {
        cloud[i].x = (double)(rand())/RAND_MAX*(max - min) + min;
        cloud[i].y = (double)(rand())/RAND_MAX*(max - min) + min;
    }
    return cloud;
}

void QuickSort(const vector<point>& _cloud, vector<int>& list_point, int size)
{
    //Указатели в начало и в конец
    int i = 1;
    int j = size - 1;

    //  Выбранный элемент
    point mid = _cloud[size / 2];

    //  Делим массив
    do {
        //  В левой части оставляем, меньше выбранного
        while (SideSpace(_cloud[list_point[0]], _cloud[list_point[i]], mid) < 0)
            i++;

        //  В правой части оставляем, больше выбранного
        while (SideSpace(_cloud[list_point[0]], _cloud[list_point[j]], mid) > 0)
            j--;

        //  Меняем элементы местами
        if (i <= j) {
            swap(list_point[i], list_point[j]);

            i++;
            j--;
        }
    } while (i <= j);


    //  Рекурсивные вызовы
    if (j > 0)
        QuickSort(_cloud, list_point, j + 1);

    if (i < size)
        QuickSort(_cloud, list_point, size - i);
}

vector<int> SequentialPassageGraham(vector<point>& cloud)
{
    if (cloud.size() == 1) {
        vector<int> res;
        res.push_back(0);
        return res;
    }

    //  Первый шаг - найдем точку точно входящую в МВО
    //  Точка с минимальной левой координатой х, если несколько, то и у

    vector<int> number_point(cloud.size());
    for (int i = 0; i < cloud.size(); i++)
        number_point[i] = i;

    for (int i = 1; i < cloud.size(); i++)
        if (cloud[number_point[i]].x < cloud[number_point[0]].x)
            swap(number_point[0], number_point[i]);

    //  Второй шаг - сортировка точек

    //QuickSort(cloud, number_point, cloud.size());
    InsertionSort(cloud, number_point);

    //  Третий шаг - убрать ненужные вершины( правый поворот плохой)

    //  2 элеманта точно войдут в оболочку
    vector<int> mch;
    mch.push_back(number_point[0]);
    mch.push_back(number_point[1]);

    //  если правый поворот, то срезаем вершину
    //  левый, добавляем
    for (int i = 2; i < cloud.size(); i++) {
        while (SideSpace(cloud[mch[mch.size() - 2]], cloud[mch[mch.size() - 1]], cloud[number_point[i]]) < 0) {
            mch.pop_back();
        }
        mch.push_back(number_point[i]);
    }
    // номера точек в исходном массиве(облаке)
    return mch;
}

vector<int> ParallelPassageGraham(vector<point>& cloud) {

    int RANK;
    int SIZE;
    MPI_Comm_size(MPI_COMM_WORLD, &SIZE);
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);

    if (SIZE == 1) {
        return SequentialPassageGraham(cloud);
    }
    //  Перевили облако в массив даблов
    vector<double> cloud_double;
    for (int i = 0; i < cloud.size(); i++) {
        cloud_double.push_back(cloud[i].x);
        cloud_double.push_back(cloud[i].y);
    }

    int div = cloud.size() / (SIZE - 1);
    int mod = cloud.size() % (SIZE - 1);
    vector<double> part_cloud(2 * div);

    if (RANK == 0) {
        for (int proc = 0; proc < SIZE - 1; proc++) {
            if (proc + 1 <= mod) {
                MPI_Send(&cloud_double[0] + (div + 1) * 2 * proc, 2 * (div + 1), MPI_DOUBLE, proc + 1, 0, MPI_COMM_WORLD);
            }
            else {
                MPI_Send(&cloud_double[0] + 2 * (div * proc + mod), 2 * div, MPI_DOUBLE, proc + 1, 0, MPI_COMM_WORLD);
            }
        }
    }
    else {
        MPI_Status stat;

        if (RANK <= mod) {
            part_cloud.resize(2 * (div + 1));
            MPI_Recv(&part_cloud[0], part_cloud.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);
        }
        else {
            MPI_Recv(&part_cloud[0], part_cloud.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);
        }


        // Перевели частицу облака из даблов в точки
        vector<point> cloudlet(part_cloud.size() / 2);

        for (int i = 0; i < cloudlet.size(); i++) {
            cloudlet[i].x = part_cloud[2 * i];
            cloudlet[i].y = part_cloud[2 * i + 1];
        }

        //  Это оболочка в списке точек маленького облачка
        vector<int> partical_shell_point = SequentialPassageGraham(cloudlet);

        // А здесь переведем нумерацию оболочки, как в основном облаке
        for (int i = 0; i < partical_shell_point.size(); i++) {
            for (int j = 0; j < cloud.size(); j++) {
                if (cloudlet[partical_shell_point[i]].x == cloud[j].x && cloudlet[partical_shell_point[i]].y == cloud[j].y) {
                    partical_shell_point[i] = j;
                    break;
                }
            }
        }

        MPI_Send(&partical_shell_point[0], partical_shell_point.size(), MPI_INT, 0, 0, MPI_COMM_WORLD);

    }

    //  Куда принять номера точек оболочки
    vector<int> shells;

    //  Все оболочки вместе, для передачи в последовательную функцию
    vector<point> big_shell;

    if (RANK == 0) {

        for (int i = 1; i < SIZE; i++) {
            MPI_Status stat;
            ////  Приняли длину оболочки
            int NumElems;
            MPI_Probe(i, 0, MPI_COMM_WORLD, &stat);
            MPI_Get_count(&stat, MPI_INT, &NumElems);
            shells.resize(NumElems);

            // Прияняли саму оболочку в виде номеров точек в большом облаке
            MPI_Recv(&shells[0], NumElems, MPI_INT, i, 0, MPI_COMM_WORLD, &stat);

            //  Перевели номера точек в сами точки в большом списке
            for (int j = 0; j < shells.size(); j++) {
                big_shell.push_back(cloud[shells[j]]);
            }
        }

        //  Из большой и повторяющейся оболочки сделали нормальную( вернулась в виде номеров точек)
        

        shells = SequentialPassageGraham(big_shell);

        //  Восстановили нумерацию, как в большом облаке и вернули этот список номеров точек
        for (int i = 0; i < shells.size(); i++) {
            for (int j = 0; j < cloud.size(); j++) {
                if (big_shell[shells[i]].x == cloud[j].x && big_shell[shells[i]].y == cloud[j].y) {
                    shells[i] = j;
                    break;
                }
            }
        }


        return shells;
    }


}

int main(int argc, char** argv)
{
    int RANK;

    vector<point> cloud = getRandomCloud(10000);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
    vector<int> list_parall = ParallelPassageGraham(cloud);

    if (RANK == 0) {
        //cout << "Parallel  ---  ";
        //for (int i = 0; i < list_parall.size(); i++)
        //{
        //    cout << list_parall[i] << "  ";
        //}
        //cout << endl;

        //vector<int> list_sequen = SequentialPassageGraham(cloud);
        //cout << "Sequential  ---  ";
        //for (int i = 0; i < list_sequen.size(); i++)
        //{
        //    cout << list_sequen[i] << "  ";
        //}

        vector<int> list_sequen = SequentialPassageGraham(cloud);

        if (list_parall == list_sequen) {
            cout << "TRUE" << endl;
        }
        else {
            cout << "FALSE" << endl;
        }
    }
    MPI_Finalize();

    return 0;
}


