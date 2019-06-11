// Copyright 2019 Lobanov Andrey

#define _SCL_SECURE_NO_WARNINGS
#define namount 1048576

#include <tbb/tbb.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

void GenerateArr(int *const arr, int const len) {
    srand((unsigned int)time(NULL));

    for (int i = 0; i < len; i++) {
        arr[i] = std::rand() % 1000;
    }
}

void ShowArray(const int* arr, const int len) {
    if (len > 65) {
        std::cout << "size is too big to display" << std::endl;
        return;
    }
    for (int i = 0; i < len; i++)
        std::cout << arr[i] << " ";

    std::cout << std::endl;
}

bool SortCheck(int const *arr, int len) {
    for (int i = 0; i < len - 1; i++)
        if (arr[i] > arr[i + 1])
            return false;
    return true;
}

const int GetMax(int *const arr, int const len) {
    int max = arr[0];
    for (int i = 1; i < len; i++) {
        if (arr[i] > max) {
             max = arr[i];
        }
    }
    return max;
}

void CountSort(int *const arr, int len, int exp) {
    int* output = new int[len];
    int i, count[10] = { 0 };

    for (i = 0; i < len; i++) {
        count[(arr[i] / exp) % 10]++;
    }
    for (i = 1; i < 10; i++) {
        count[i] += count[i - 1];
    }
    for (i = len - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
    }

    for (i = 0; i < len; i++) {
        arr[i] = output[i];
    }
}

void RadixSort(int *const arr, int len) {
    int exp, m;
    m = GetMax(arr, len);

    for (exp = 1; m / exp > 0; exp *= 10) {
        CountSort(arr, len, exp);
    }
}

void BatcherMerge(int *arr, int l, int r) {
    int count = r - l + 1;

    for (int k = count / 2; k > 0; k /= 2)
        for (int j = k % (count / 2); j + k < count; j += k + k)
            for (int i = 0; i < k; i++)
                if (arr[l + j + i] > arr[l + j + i + k]) {
                    int temp = arr[l + j + i];
                    arr[l + j + i] = arr[l + j + i + k];
                    arr[l + j + i + k] = temp;
                }
}

void RadixWBatchermergeTbb(int *arr, int len, int amount, int threads) {
    int step = 1;
    int exp = 1;

    tbb::parallel_for(tbb::blocked_range<int>(0, threads, 1), [=](tbb::blocked_range<int>& r) {
        int begin = r.begin();
        int end = r.end();

    for (int i = begin; i != end; ++i)
            RadixSort(arr + i * amount, amount);
    });

    while (step < threads) {
        step = static_cast<int>(pow(2, exp));

        tbb::parallel_for(tbb::blocked_range<int>(0, threads / step, 1), [=](tbb::blocked_range<int>& r) {
            int begin = r.begin();
            int end = r.end();

            for (int i = begin; i != end; i++)
                BatcherMerge(arr, i * amount * step, (step * i * amount + step * amount) - 1);
        });
        exp++;
    }
}

int main(int argc, char* argv[]) {
    int len = namount;

    double dt, dt2;

    int threads = tbb::task_scheduler_init::default_num_threads();
    tbb::task_scheduler_init init(threads);

    int *arr = new int[len];
    int *arr2 = new int[len];

    GenerateArr(arr, len);
    GenerateArr(arr2, len);

    int shift = len / threads;

    tbb::tick_count t1 = tbb::tick_count::now();
    RadixWBatchermergeTbb(arr, len, shift, threads);
    tbb::tick_count t2 = tbb::tick_count::now();

    dt = (t2 - t1).seconds();

    t1 = tbb::tick_count::now();
    RadixSort(arr2, len);
    t2 = tbb::tick_count::now();

    dt2 = (t2 - t1).seconds();

    ShowArray(arr, len);

    if (SortCheck(arr, len))
        std::cout << "\nArray is sorted";
    else
        std::cout << "\nUnsorted!!!";

    std::cout << "\n\nElements: " << namount;
    std::cout << "\n\nTime PP for [" << threads << "] threads : " << dt;
    std::cout << "\n\nTime N_PP :" << dt2;
    std::cout << "\n\nBoost: " << dt2 / dt << std::endl;

    delete[] arr;
    delete[] arr2;

    return 0;
}
