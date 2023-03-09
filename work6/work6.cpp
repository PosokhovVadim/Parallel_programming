#include <omp.h>
#include <iostream>
#include <time.h>
#include <cmath>
#include <time.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

void PrintMatrix (double** pMatrix, int RowCount, int ColCount) {
    int i, j; // Loop variables
    for (i=0; i<RowCount; i++) {
        for (j=0; j<ColCount; j++)
            cout << pMatrix[i][j];
        printf("\n");
    }
}

void PrintVector (double pVector[], int Size) {
    int i;
    for (i=0; i<Size; i++)
        printf("%7.4f ", pVector[i]);
}

// Function that gets the timestamp in seconds

double CountRow(int from, int i, int n, double** matrix, double* res_vector){
    double res_of_sum;

    for (int j = from; j < n; j++){
        res_of_sum += matrix[i][j] * res_vector[j];
    }

    return res_of_sum;

}
double* RelaxationMethod(double** matrix, double b_vector[], double res_vector[], int n, double omega) {
	int i, j, k = 0;
    for (int i = 0; i < n; i++) {
        res_vector[i] = 0;
    }
    for (int i = 0; i < n; i++) {
            double sum = 0.0;
            // for (int j = 0; j < n; j++) {
            //     if (i != j) {
            //         sum += matrix[i][j] * res_vector[j];
            //     }
            // }

            double new_x = (1 - omega) * res_vector[i] + (omega * ( (b_vector[i] - CountRow(1,i, i - 1, matrix, res_vector) 
                                                                               - CountRow(i + 1,i ,n, matrix, res_vector)) / matrix[i][i] )); 
            res_vector[i] = new_x;
        }
    return res_vector;
}
int main()
{
	setlocale(LC_CTYPE, "RUSSIAN");
    int n;
    int omega;
    cout << "Введите размерность матрицы:";
    cin >> n;
    cout << endl;
    cout << "Введите параметр релаксации: ";
    cin >> omega;
    cout << endl;

    double** matrix = new double* [n];
    for (int i = 0; i < n; i++)
        matrix[i] = new double [n];
    double* vector = new double[n];

    double* b_vector = new double[n];
    //Заполнение
    double x;
    for (int i =0; i < n; i++) {
        for (int j =0; j < n; j++) {
            cin >> x;
            matrix[i][j] = x;
        }
    }
    for (int i = 0; i < n; i++){
        cin >> x;
        b_vector[i] = x; 
    } 
    

    printf("Initial Matrix \n");
    PrintMatrix(matrix, n, n);
    
    steady_clock::time_point start = steady_clock::now();
    RelaxationMethod(matrix, vector,b_vector, n, omega);
    steady_clock::time_point finish = steady_clock::now();
    duration<double> diff = finish - start;
    
    
    //The matrix and the vector output
    printf("\nResult Vector \n");
    PrintVector(vector, n);
    cout <<"\n Time of execution: " << diff.count() << endl;
    //printf("\n Time of execution: %lld\n", diff.count());
 // Computational process termination
    return 0;


}