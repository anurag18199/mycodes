#include <iostream>
#include <omp.h>
#include <vector>
#include <limits>

using namespace std;

int main() {
    const int N = 10000;
    vector<int> data(N);

    // Initialize the data
    for (int i = 0; i < N; ++i)
        data[i] = i % 1000;  // example data

    int min_val = numeric_limits<int>::max();
    int max_val = numeric_limits<int>::min();
    long long sum = 0;

    #pragma omp parallel for reduction(min:min_val) reduction(max:max_val) reduction(+:sum)
    for (int i = 0; i < N; ++i) {
        int val = data[i];
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
        sum += val;
        // cout<<val<<" ";
    }
    // cout<<endl;

    double avg = (double)sum / N;

    cout << "Min: " << min_val << endl;
    cout << "Max: " << max_val << endl;
    cout << "Sum: " << sum << endl;
    cout << "Average: " << avg << endl;

    return 0;
}
