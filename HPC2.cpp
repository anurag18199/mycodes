#include <iostream>
#include <vector>
#include <random>
#include <omp.h>

using namespace std;

// Sequential Bubble Sort
void bubbleSort(vector<int>& a) {
    int n = a.size();
    for (int i = 0; i < n - 1; ++i)
        for (int j = 0; j < n - i - 1; ++j)
            if (a[j] > a[j+1])
                swap(a[j], a[j+1]);
}

// Parallel Bubble Sort via Odd–Even Transposition
void parallelBubbleSort(vector<int>& a) {
    int n = a.size();
    for (int phase = 0; phase < n; ++phase) {
        int start = phase % 2; // even or odd phase

        #pragma omp parallel for default(shared) schedule(static)
        for (int j = start; j < n - 1; j += 2) {
            if (a[j] > a[j+1])
                swap(a[j], a[j+1]);
        }
    }
}


// Merge routine for Merge Sort
void merge(vector<int>& a, int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;
    vector<int> L(n1), R(n2);
    for (int i = 0; i < n1; ++i) L[i] = a[l + i];
    for (int j = 0; j < n2; ++j) R[j] = a[m + 1 + j];
    int i = 0, j = 0, k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) a[k++] = L[i++];
        else             a[k++] = R[j++];
    }
    while (i < n1) a[k++] = L[i++];
    while (j < n2) a[k++] = R[j++];
}

// Sequential Merge Sort
void mergeSort(vector<int>& a, int l, int r) {
    if (l >= r) return;
    int m = l + (r - l) / 2;
    mergeSort(a, l, m);
    mergeSort(a, m+1, r);
    merge(a, l, m, r);
}

// Parallel Merge Sort using OpenMP tasks
void parallelMergeSort(vector<int>& a, int l, int r, int depth = 0) {
    const int THRESHOLD = 1000; 
    if (r - l < THRESHOLD) {
        // small enough — do serial
        mergeSort(a, l, r);
        return;
    }
    int m = l + (r - l) / 2;
    #pragma omp task default(shared) firstprivate(l, m, depth)
    parallelMergeSort(a, l, m, depth+1);
    #pragma omp task default(shared) firstprivate(m, r, depth)
    parallelMergeSort(a, m+1, r, depth+1);
    #pragma omp taskwait
    merge(a, l, m, r);
}

int main() {
    const int N = 10000;                // array size
    vector<int> original(N);

    // initialize RNG
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(0, 1000000);
    for (int i = 0; i < N; ++i)
        original[i] = dist(gen);

    double t1, t2;

    // 1) Sequential Bubble Sort
    vector<int> a1 = original;
    t1 = omp_get_wtime();
    bubbleSort(a1);
    t2 = omp_get_wtime();
    cout << "Sequential Bubble Sort: " << (t2 - t1) << " s\n";

    // 2) Parallel Bubble Sort
    vector<int> a2 = original;
    t1 = omp_get_wtime();
    parallelBubbleSort(a2);
    t2 = omp_get_wtime();
    cout << "Parallel Bubble Sort:   " << (t2 - t1) << " s\n";

    // 3) Sequential Merge Sort
    vector<int> a3 = original;
    t1 = omp_get_wtime();
    mergeSort(a3, 0, N-1);
    t2 = omp_get_wtime();
    cout << "Sequential Merge Sort:  " << (t2 - t1) << " s\n";

    // 4) Parallel Merge Sort
    vector<int> a4 = original;
    t1 = omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp single
        parallelMergeSort(a4, 0, N-1);
    }
    t2 = omp_get_wtime();
    cout << "Parallel Merge Sort:    " << (t2 - t1) << " s\n";

    return 0;
}
