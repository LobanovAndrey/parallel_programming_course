//  Copyright 2019 Lobanov Andrey
#define nSize 1000
#define namount 1048576

#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

void generateArr(int *const arr, int const len) {
    srand((unsigned int)time(NULL));

    for (int i = 0; i < len; i++) {
        arr[i] = std::rand() % nSize;
    }
}

void showArray(const int* arr, const int length) {
    if (length > 65) {
        std::cout << "size is too big to display" << std::endl;
        return;
    }
    for (int i = 0; i < length; i++)
        std::cout << arr[i] << " ";

    std::cout << std::endl;
}

bool sortCheck(int const *arr, int len) {
    for (int i = 0; i < len - 1; i++)
        if (arr[i] > arr[i + 1])
           return false;
    return true;
}

const int getMax(int *const arr, int const len) {
    int max = arr[0];
    for (int i = 1; i < len; i++) {
        if (arr[i] > max) {
            max = arr[i];
        }
    }
    return max;
}

void countSort(int *const arr, int len, int exp) {
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

void radixSort(int *const arr, int len) {
    int exp, m;
    m = getMax(arr, len);

    for (exp = 1; m / exp > 0; exp *= 10) {
        countSort(arr, len, exp);
    }
}

void merge(int *arr, int l, int r) {
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

int main(int argc, char* argv[]) {
    int threads = 4;

    double t1, t2, dt, dt2;
    int len = namount;
    int *arr = new int[len];
    int *arr2 = new int[len];

    int step = 1;
    int exp = 1;

    generateArr(arr, len);
    generateArr(arr2, len);

    int* shift = new int[threads];
    int* amount = new int[threads];

    t1 = omp_get_wtime();

#pragma omp parallel shared(arr, shift, amount, exp) num_threads(threads)
    {
        int t_id = omp_get_thread_num();

        shift[t_id] = t_id * (len / threads);

        if (t_id == threads - 1)
            amount[t_id] = (len / threads) + (len % threads);
        else
            amount[t_id] = len / threads;

        radixSort(arr + shift[t_id], amount[t_id]);

        while (step < threads) {
#pragma omp barrier
            step = static_cast<int>(pow(2, exp));

            if (t_id < (threads / step)) {
                merge(arr, shift[t_id] * step, (step * shift[t_id] + step * amount[t_id]) - 1);
            }
#pragma omp barrier
#pragma omp single
            {
                exp++;
            }
        }
    }

    t2 = omp_get_wtime();
    dt = t2 - t1;

    t1 = omp_get_wtime();
    radixSort(arr2, len);
    t2 = omp_get_wtime();
    dt2 = t2 - t1;

    showArray(arr, len);

    if (sortCheck(arr, len))
        std::cout << "\nSorted";
    else
        std::cout << "\nUnsorted!!!";

    std::cout << "\n\nTime P :" << dt;
    std::cout << "\n\nTime NP :" << dt2;

    delete[] arr;
    delete[] arr2;
    delete[] shift;
    delete[] amount;

    return 0;
}
