#include <iostream>
#include <omp.h>
#include <cmath>
#include  <time.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

double func(double x){
    return 1.0 / (1.0 + x * x);
}

double Integral_Rectangle_Method(const double a, const double b, const double h){
    int n = (int)((b - a)/ h);
    double sum = 0.0;
    int i;
    double x;
    #pragma omp parallel for private(x) reduction(+:sum)
    for(i = 0; i < n; i++){
        x = a + i*h + h/2;
        sum += func(x);
    }
    return h * sum;
    
}


double Integral_Simpson_Method(const double a, const double b, const double h){
    int n = (int)((b - a)/ (2.0 * h));
    int i;
    double x;
    double f_a = func(a);
    double f_b = func(b);
    double sum_1 = 0.0;
    double sum_2 = 0.0;
    #pragma omp parallel for private(x) reduction(+:sum_1)
    for (int i = 1; i < n; i += 2)
    {
        x = a + i*h;
        sum_1 += func(x);
    }
    #pragma omp parallel for private(x) reduction(+:sum_2)
    for (int i = 2; i < n; i += 2)
    {
        x = a + i*h;
        sum_2 += func(x);
    }

    return (f_a + f_b + 4.0 * sum_1 + 2.0 * sum_2) * h / 3.0;
}


double experiment(double *res, double f(double, double, double)){
    double stime, ftime; // время начала и конца расчета
    double a, b, h;
    a = 0.0;
    b = 100000;
    h = 0.001; 
    steady_clock::time_point start = steady_clock::now();
    *res = f(a, b, h);

    steady_clock::time_point end = steady_clock::now();
    duration<double> diff = end - start;
    // cout << "\nTime = " << diff.count() << endl;

    // stime = clock( );
    // ftime = clock( );
    return diff.count();
}


int main() {

    double res;
    int numbExp = 10;
    double time;
    double min_time; 
    double max_time; 
    double avg_time;
    double(*func_pointer)(double a, double b, double h);
    int k;
    while (true)
    {
        cout << "0 - Rectangle_Method\n1 - Simpson_Method\n";
        cin >> k;
        switch (k){
            case 0: func_pointer = Integral_Rectangle_Method; break;
            case 1: func_pointer = Integral_Simpson_Method; break;
            default: return 0;
        }
        min_time = max_time = avg_time = experiment(&res, func_pointer);
        for(int i = 0; i < numbExp - 1; i ++) {
            time = experiment(&res, func_pointer);
            avg_time += time;
            if(max_time < time) max_time = time;
            if(min_time > time) min_time = time;
        }
        cout << "execution time : " << avg_time / numbExp << "; " <<
        min_time<< "; " << max_time << endl;
        cout.precision(8);
        cout << "integral value : " << res << endl; 
    }
    return 0;

}
