// Copyright 2020 Bezrodnov Dmitry
#ifndef MODULES_TASK_3_BEZRODNOV_D_PASSAGE_GRAHAM_PASSAGE_GRAHAM_H_
#define MODULES_TASK_3_BEZRODNOV_D_PASSAGE_GRAHAM_PASSAGE_GRAHAM_H_

#include <vector>

struct point {
    double x;
    double y;

    point() {
        x = 0;
        y = 0;
    }

    //  Конструктор копироваия
    point(const point& op2) {
        x = op2.x;
        y = op2.y;
    }

    point(double _x, double _y) {
        x = _x;
        y = _y;
    }

    //  Перегрузка присваивания
    point& operator=(const point& op2) {
        x = op2.x;
        y = op2.y;
        return *this;
    }
};

double SideSpace(point A1, point A2, point B);

std::vector<point> getRandomCloud(int size);

point getRandomPoint();

std::vector<int> InsertionSort(const std::vector<point>& _cloud, const std::vector<int>& list_point);

std::vector<int> SequentialPassageGraham(const std::vector<point>& cloud);

std::vector<int> ParallelPassageGraham(const std::vector<point>& cloud);

#endif  //  MODULES_TASK_3_BEZRODNOV_D_PASSAGE_GRAHAM_PASSAGE_GRAHAM_H_

