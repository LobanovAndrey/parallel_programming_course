//  Copyright 2019 Lobanov Andrey
#define nSize 1000
#define nAmount 4096

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

void ShowArray(const int* arr, const int length) {
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

void compexch(int &const a, int &const b)
{
    if (b < a)
        std::swap(a, b);
}

void radixSort(int *const arr, int len) {
    int exp, m;
    m = getMax(arr, len);

    for (exp = 1; m / exp > 0; exp *= 10) {
        countSort(arr, len, exp);
    }
}

void merge(int *const arr, int l, int r)
{
    int count = r - l + 1;

    for (int k = count / 2; k > 0; k /= 2)
        for (int j = k % (count / 2); j + k < count; j += k + k)
            for (int i = 0; i < k; i++)
                compexch(arr[l + j + i], arr[l + j + i + k]);
}

int main(int argc, char* argv[])
{
    int threads = 4;
    int t_id, shift, amount;
    double t1, t2, dt, dt2;
    int len = nAmount;
    int	*arr = new int[len];
    int *arr2 = new int[len];

    int step = 1;
    int exp = 1;

    generateArr(arr, len);
    generateArr(arr2, len);

    t1 = omp_get_wtime();

#pragma omp parallel shared(arr, exp) private(t_id, shift, amount) num_threads(threads) 
    {
        t_id = omp_get_thread_num();

        shift = t_id * (len / threads);

        if (t_id == threads - 1)
            amount = (len / threads) + (len % threads);
        else
            amount = len / threads;

        radixSort(arr + shift, amount);

        while (step < threads) {
#pragma omp barrier
            step = pow(2, exp);

            if (t_id < (threads / step)) {
                merge(arr, shift * step, (step * shift + step * amount) - 1);
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

    ShowArray(arr, len);

    if (sortCheck(arr, len))
        std::cout << "\nSorted";
    else
        std::cout << "\nUnsorted!!!";

    std::cout << "\nTime P :" << dt;
    std::cout << "\nTime NP :" << dt2;

    delete[] arr;
    delete[] arr2;

    return 0;
}
