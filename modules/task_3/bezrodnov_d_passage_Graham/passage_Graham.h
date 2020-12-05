// Copyright 2020 Bezrodnov Dmitry
#ifndef MODULES_TASK_3_BEZRODNOV_D_PASSAGE_GRAHAM_PASSAGE_GRAHAM_H_
#define MODULES_TASK_3_BEZRODNOV_D_PASSAGE_GRAHAM_PASSAGE_GRAHAM_H_
#include <vector>

struct point {
    double x;
    double y;
};

double SideSpace(point A1, point A2, point B);

std::vector<point> getRandomCloud(int size);

void InsertionSort(const std::vector<point>& _cloud, std::vector<int>& list_point);

std::vector<int> SequentialPassageGraham(const std::vector<point>& cloud);

std::vector<int> ParallelPassageGraham(const std::vector<point>& cloud);

#endif  //  MODULES_TASK_3_BEZRODNOV_D_PASSAGE_GRAHAM_PASSAGE_GRAHAM_H_

