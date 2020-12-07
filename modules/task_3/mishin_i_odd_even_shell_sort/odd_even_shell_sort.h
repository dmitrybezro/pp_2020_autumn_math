// Copyright 2020 Mishin Ilya
#ifndef MODULES_TASK_3_MISHIN_I_ODD_EVEN_SHELL_SORT_ODD_EVEN_SHELL_SORT_H_
#define MODULES_TASK_3_MISHIN_I_ODD_EVEN_SHELL_SORT_ODD_EVEN_SHELL_SORT_H_

#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <utility>

template< typename RandomAccessIterator, typename Compare >
void shell_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
    for (auto d = (last - first) / 2; d != 0; d /= 2)
        for (auto i = first + d; i != last; ++i)
            for (auto j = i; j - first >= d && comp(*j, *(j - d)); j -= d)
                std::swap(*j, *(j - d));
}

std::vector<int> getRandomVector(int size);

int* parallel_shell_sort(int* container, int array_size);

#endif  // MODULES_TASK_3_MISHIN_I_ODD_EVEN_SHELL_SORT_ODD_EVEN_SHELL_SORT_H_
