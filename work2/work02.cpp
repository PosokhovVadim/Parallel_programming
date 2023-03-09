#include <iostream>
#include <omp.h>
#include <cmath>
#include  <time.h>
#include <chrono>

using namespace std;
using namespace std::chrono;
#define EXP 2.718281828459045235360287471352
#define PI 3.1415926535897932384626433832795

double func(double x, double y){
    return pow(EXP, (sin(PI * x) * cos(PI * y))) + 1.0;
}

double Integral_Rectangle_Method(const double a1, const double b1, const double h1, const double a2, const double b2, const double h2){
    int n = (int)((b1 - a1)/ h1);
    int m = (int)((b2-a2)/ h2);
    double sum = 0.0;
    int i, j;
    double x,y;
    #pragma omp parallel for private(x) reduction(+ : sum)
    for(i = 0; i < n; i++){
        x = a1 + i*h1 + h1/2.0;
        #pragma omp parallel for private(y) 
        for(j = 0; j < m; j++){
            y = a2 + j*h2 + h2/2.0;
            sum += h1 * h2 *  func(x,y) /  ((b1 - a1) * (b2 -a2));
        }
    }
    return sum;
    
}


double experiment(double *res){
    double stime, ftime; // время начала и конца расчета
    //omp_set_num_threads(8);
    double a1, b1, h1, a2, b2, h2;
    a1 = a2 = 0.0;
    b1 = b2 = 16.0;
    h1 = h2 = 0.001; 
    // steady_clock::time_point start = steady_clock::now();
    //*res = Integral_Rectangle_Method(a1,b1,h1, a2, b2, h2);

    // steady_clock::time_point end = steady_clock::now();
    // duration<double> diff = end - start;
    // cout << "\nTime = " << diff.count() << endl;

    stime = clock( );
    *res = Integral_Rectangle_Method(a1,b1,h1, a2, b2, h2);
    ftime = clock( );
    return (ftime - stime) / CLOCKS_PER_SEC;
}


int main() {
    double res;
    int numbExp = 10;
    double time;
    double min_time; 
    double max_time; 
    double avg_time;
    //double(*func_pointer)(double a, double b, double h);
    min_time = max_time = avg_time = experiment(&res);
    for(int i = 0; i < numbExp - 1; i ++) {
        time = experiment(&res);
        avg_time += time;
        if(max_time < time) max_time = time;
        if(min_time > time) min_time = time;
    }
    cout << "execution time : " << avg_time / numbExp << "; " <<
    min_time<< "; " << max_time << endl;
    cout.precision(8);
    cout << "integral value : " << res << endl; 
        
    return 0;

}
