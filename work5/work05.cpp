#include <iostream>
#include <omp.h>
#include <cmath>
#include  <time.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

int* pPivotPos;
int* pPivotIter;

typedef struct {
    int PivotRow;
    double MaxValue;
} 
TThreadPivotRow;


void PrintMatrix (double* pMatrix, int RowCount, int ColCount) {
    int i, j; // Loop variables
    for (i=0; i<RowCount; i++) {
        for (j=0; j<ColCount; j++)
            printf("%7.4f ", pMatrix[i*RowCount+j]);
        printf("\n");
    }
}


void PrintVector (double* pVector, int Size) {
    int i;
    for (i=0; i<Size; i++)
        printf("%7.4f ", pVector[i]);
}



void DataInitialization (double* pMatrix, double* pVector, int Size) {
    int i, j; // Loop variables
    //srand(unsigned(clock()));
    cout << "Коэффициенты матрицы: "<< endl;
    for (i=0; i<Size; i++) {
        for (j=0; j<Size; j++) {
           cin >> pMatrix[i*Size+j];
        }
        cout << "\n";
    }
    cout << "Коэффициенты результата:" << endl;
    for (int i = 0; i < Size; i++ ){
        cin >> pVector[i];
    }

}



// Back substation
void ParallelBackSubstitution (double* pMatrix, double* pVector, double* pResult, int Size) {
    int RowIndex, Row;
    for (int i=Size-1; i>=0; i--) {
        RowIndex = pPivotPos[i];
        pResult[i] = pVector[RowIndex]/pMatrix[Size*RowIndex+i];
    #pragma omp parallel for private (Row)
    for (int j=0; j<i; j++) {
        Row = pPivotPos[j];
        pVector[Row] -= pMatrix[Row*Size+i]*pResult[i];
        pMatrix[Row*Size+i] = 0;
        }
    }
}


void ProcessInitialization (double* &pMatrix, double* &pVector, double* &pResult, int &Size) {
 // Setting the size of the matrix and the vector
    do {
        printf("\nEnter size of the matrix and the vector: ");
        scanf("%d", &Size);
        printf("\nChosen size = %d \n", Size);
    if (Size <= 0)
        printf("\nSize of objects must be greater than 0!\n");
    } while (Size <= 0);
 // Memory allocation
    pMatrix = new double [Size*Size];
    pVector = new double [Size];
    pResult = new double [Size];
 // Initialization of the matrix and the vector elements

    //DummyDataInitialization(pMatrix, pVector, Size);
    DataInitialization(pMatrix, pVector, Size);
}

// Column elimination
void ParallelColumnElimination (double* pMatrix, double* pVector, int Pivot, int Iter, int Size) {
    double PivotValue, PivotFactor;
    PivotValue = pMatrix[Pivot*Size+Iter];
    #pragma omp parallel for private(PivotFactor) schedule(dynamic,1)
    for (int i=0; i<Size; i++) {
        if (pPivotIter[i] == -1) {
            PivotFactor = pMatrix[i*Size+Iter] / PivotValue;
            for (int j=Iter; j<Size; j++) {
                pMatrix[i*Size + j] -= PivotFactor * pMatrix[Pivot*Size+j];
            }
        pVector[i] -= PivotFactor * pVector[Pivot];
        }
    }
}


int ParallelFindPivotRow(double* pMatrix, int Size, int Iter) {
    int PivotRow = -1; // The index of the pivot row
    int MaxValue = 0; // The value of the pivot element
    int i; // Loop variable
    #pragma omp parallel
    {
    TThreadPivotRow ThreadPivotRow;
    ThreadPivotRow.MaxValue = 0;
    ThreadPivotRow.PivotRow = -1;
    #pragma omp for
    for (i=0; i<Size; i++) {
        if ((pPivotIter[i] == -1) && (fabs(pMatrix[i*Size+Iter]) > ThreadPivotRow.MaxValue)) {
            ThreadPivotRow.PivotRow = i;
            ThreadPivotRow.MaxValue = fabs(pMatrix[i*Size+Iter]);
        }
    }
    //printf("\n Local thread (id = %i) pivot row : %i " ,omp_get_thread_num(), ThreadPivotRow.PivotRow);
    
    #pragma omp critical
    {
    if (ThreadPivotRow.MaxValue > MaxValue){
        MaxValue = ThreadPivotRow.MaxValue;
        PivotRow = ThreadPivotRow.PivotRow;
    }
    
    } // pragma omp critical// pragma omp parallel
    }
    return PivotRow;
}

void ParallelGaussianElimination(double* pMatrix,double* pVector, int Size) {
    int Iter; // The number of the iteration of the Gaussian
    // elimination
    int PivotRow; // The number of the current pivot row
    for (Iter=0; Iter<Size; Iter++) {
    // Finding the pivot row
        PivotRow = ParallelFindPivotRow(pMatrix, Size, Iter);
        pPivotPos[Iter] = PivotRow;
        pPivotIter[PivotRow] = Iter;
        ParallelColumnElimination(pMatrix, pVector, PivotRow, Iter, Size);
    }
}


void ParallelResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size) {
// Memory allocation
    pPivotPos = new int [Size];
    pPivotIter = new int [Size];
    for (int i=0; i<Size; i++) {
        pPivotIter[i] = -1;
    }
    ParallelGaussianElimination(pMatrix, pVector, Size);
    ParallelBackSubstitution (pMatrix, pVector, pResult, Size);
    delete [] pPivotPos;
    delete [] pPivotIter;
 }


void ProcessTermination (double* pMatrix, double* pVector, double* pResult) {
    delete [] pMatrix;
    delete [] pVector;
    delete [] pResult;
} 


int main(){

    double* pMatrix; // The matrix of the linear system
    double* pVector; // The right parts of the linear system
    double* pResult; // The result vector
    int Size; // The sizes of the initial matrix and the vector
    //printf("Serial Gauss algorithm for solving linear systems\n");
    ProcessInitialization(pMatrix, pVector, pResult, Size);
    printf("Initial Matrix \n");
    PrintMatrix(pMatrix, Size, Size);
    printf("Initial Vector \n");
    PrintVector(pVector, Size);
    
    steady_clock::time_point start = steady_clock::now();
    //start = omp_get_wtime();
    ParallelResultCalculation(pMatrix, pVector, pResult, Size);
    //finish = omp_get_wtime();
    steady_clock::time_point finish = steady_clock::now();
    duration<double> diff = finish - start;
    
    
    
    
    //The matrix and the vector output
    printf ("Recalculated Matrix \n");
    PrintMatrix(pMatrix, Size, Size);
    printf("Recalculated Vector \n");
    PrintVector(pVector, Size);
    printf("\nResult Vector \n");
    PrintVector(pResult, Size);
    cout <<"\n Time of execution: " << diff.count() << endl;
    //printf("\n Time of execution: %lld\n", diff.count());
 // Computational process termination
    ProcessTermination(pMatrix, pVector, pResult);
    return 0;

}