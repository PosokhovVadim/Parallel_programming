#include <omp.h>
#include <iostream>
#include <time.h>
#include <cmath>
#include <time.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

void PrintMatrix (double** matrix, int n) {
    int i, j; // Loop variables
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++)
            cout << matrix[i][j] << " ";
        printf("\n");
    }
}

void PrintVector (double pVector[], int Size) {
    int i;
    for (i=0; i<Size; i++)
        printf("%7.4f ", pVector[i]);
}


double* RelaxationMethod(double** matrix, double* b_vector, double* res_vector, int n, double omega) {
	int i, j, k = 0;
    for (int i = 0; i < n; i++) {
        res_vector[i] = 1;
    }
    double sum;
    #pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i < n; i++) {
            sum = 0.0;
            #pragma omp parallel for

            for(int j=0; j<i; j++)
            {
                sum += matrix[i][j]*res_vector[j];
            }
            #pragma omp parallel for

            for(int j=i+1; j<n; j++)
            {
                sum += matrix[i][j]*res_vector[j];
            }
            double new_x = (1 - omega) * res_vector[i] + (omega * ((b_vector[i] -  sum) / matrix[i][i] )); 
            res_vector[i] = new_x;

        }
    return res_vector;
}
int main()
{
	setlocale(LC_CTYPE, "RUSSIAN");
    int n;
    srand(time(NULL));
    omp_set_num_threads(4);
    double omega;
    cout << "Введите размерность матрицы:";
    cin >> n;
    cout << endl;
    cout << "Введите параметр релаксации: ";
    cin >> omega;
    cout << endl;

    double** matrix = new double* [n];
    for (int i = 0; i < n; i++)
        matrix[i] = new double [n];
            double x;
    
    for (int i =0; i < n; i++) {
        for (int j =0; j < n; j++) {
            matrix[i][j] = rand() % 100;
        }
    }
    double* b_vector = new double[n];

        //Заполнение

    for (int i = 0; i < n; i++){
        b_vector[i] = rand() % 100; 
    } 
    double* vector = new double[n];


    //cout << "Initial matrix: " << endl;
    //PrintMatrix(matrix, n);
    
    steady_clock::time_point start = steady_clock::now();
    vector = RelaxationMethod(matrix, b_vector,vector, n, omega);
    steady_clock::time_point finish = steady_clock::now();
    duration<double> diff = finish - start;
    
    
    //The matrix and the vector output
    //printf("\nResult Vector \n");
    //PrintVector(vector, n);
    cout <<"\n Time of execution: " << diff.count() << endl;
    //printf("\n Time of execution: %lld\n", diff.count());
 // Computational process termination
    return 0;


}