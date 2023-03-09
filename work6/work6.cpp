#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <omp.h>

using namespace std;
double epsilon = 0.01;

typedef struct {
	int PivotRow; double MaxValue;
} TThreadPivotRow;
int* pSerialPivotPos;  // The Number of pivot rows selected at the//iterations 
int* pSerialPivotIter; // The Iterations, at which the rows were pivots
					   // Function for simple initialization of thematrix 
					   // and the vector elements
int* pPivotPos;  // The number of pivot rows selected at the iterations
int* pPivotIter; // The iterations, at which the rows were pivots
				 // Function for the execution of Gauss algorithm

int ParallelFindPivotRow(double* pMatrix, int Size, int Iter) {
	int PivotRow = -1;    // The index of the pivot row
	double MaxValue = 0; // The value of the pivot element
	int i;                // Loop variable

#pragma omp parallel
	{
		TThreadPivotRow ThreadPivotRow;
		ThreadPivotRow.MaxValue = 0;
		ThreadPivotRow.PivotRow = -1;
#pragma omp for
		for (i = 0; i < Size; i++) {
			if ((pPivotIter[i] == -1) &&
				(fabs(pMatrix[i * Size + Iter]) > ThreadPivotRow.MaxValue)) {
				ThreadPivotRow.PivotRow = i;
				ThreadPivotRow.MaxValue = fabs(pMatrix[i * Size + Iter]);
			}
		}
		//cout << "\n Local thread (id = %i) pivot row : %i" << omp_get_thread_num() << ThreadPivotRow.PivotRow;
#pragma omp critical
		{
			if (ThreadPivotRow.MaxValue > MaxValue) {
				MaxValue = ThreadPivotRow.MaxValue;
				PivotRow = ThreadPivotRow.PivotRow;
			}
		} // pragma omp critical
	} // pragma omp parallel
	return PivotRow;
}
// Back substation
void ParallelBackSubstitution(double* pMatrix, double* pVector,
	double* pResult, int Size) {
	int RowIndex, Row;
	for (int i = Size - 1; i >= 0; i--) {
		RowIndex = pPivotPos[i];
		pResult[i] = pVector[RowIndex] / pMatrix[Size * RowIndex + i];
#pragma omp parallel for private (Row)
		for (int j = 0; j < i; j++) {
			Row = pPivotPos[j];
			pVector[Row] -= pMatrix[Row * Size + i] * pResult[i];
			pMatrix[Row * Size + i] = 0;
		}
	}
}
// Column elimination
void ParallelColumnElimination(double* pMatrix, double* pVector,
	int Pivot, int Iter, int Size) {
	double PivotValue, PivotFactor;
	PivotValue = pMatrix[Pivot * Size + Iter];
#pragma omp parallel for private(PivotFactor) schedule(dynamic,1)
	for (int i = 0; i < Size; i++) {
		if (pPivotIter[i] == -1) {
			PivotFactor = pMatrix[i * Size + Iter] / PivotValue;
			for (int j = Iter; j < Size; j++) {
				pMatrix[i * Size + j] -= PivotFactor * pMatrix[Pivot * Size + j];
			}
			pVector[i] -= PivotFactor * pVector[Pivot];
		}
	}
}
// Gaussian elimination
void ParallelGaussianElimination(double* pMatrix, double* pVector, int Size) {
	int Iter;       // The number of the iteration of the Gaussian
					// elimination
	int PivotRow;  // The number of the current pivot row
	for (Iter = 0; Iter < Size; Iter++) {
		// Finding the pivot row
		PivotRow = ParallelFindPivotRow(pMatrix, Size, Iter);
		pPivotPos[Iter] = PivotRow;
		pPivotIter[PivotRow] = Iter;
		ParallelColumnElimination(pMatrix, pVector, PivotRow, Iter, Size);
	}
}
void ParallelResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size) {
	// Memory allocation
	pPivotPos = new int[Size];
	pPivotIter = new int[Size];
	for (int i = 0; i < Size; i++) {
		pPivotIter[i] = -1;
	}
	ParallelGaussianElimination(pMatrix, pVector, Size);
	ParallelBackSubstitution(pMatrix, pVector, pResult, Size);
	int GlobalPivotPos; GlobalPivotPos = ParallelFindPivotRow(pMatrix, Size, 0);
	//cout << "\n Global pivot row : %i", GlobalPivotPos;
	ParallelFindPivotRow(pMatrix, Size, 0);
	// Gaussian elimination
	// Back substitution 
	// Memory deallocation
	delete[] pPivotPos;
	delete[] pPivotIter;
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

void RandomDataInitialization(double* pMatrix, double* pVector, int Size) {
	int i, j;
	srand(unsigned(clock()));
	for (i = 0; i < Size; i++) {
		pVector[i] = rand() / double(1000);
		for (j = 0; j < Size; j++) {
			if (j <= i)
				pMatrix[i * Size + j] = rand() / double(1000);
			else
				pMatrix[i * Size + j] = 0;
		}
	}
}
void PRIMER(double* pMatrix, double* pVector, int Size) {
	pVector[0] = 21;
	pVector[1] = 12;
	pVector[2] = 29;
	pVector[3] = 130;
	pVector[4] = -13;

	pMatrix[0] = 5;
	pMatrix[1] = 2;
	pMatrix[2] = -7;
	pMatrix[3] = 14;
	pMatrix[4] = 0;
	pMatrix[5] = 5; 
	pMatrix[6] = -1;
	pMatrix[7] = 8;
	pMatrix[8] = -13;
	pMatrix[9] = 3;
	pMatrix[10] = 10;
	pMatrix[11] = 1;
	pMatrix[12] = -2;
	pMatrix[13] = 7;
	pMatrix[14] = -1;
	pMatrix[15] = 15;
	pMatrix[16] = 3;
	pMatrix[17] = 15;
	pMatrix[18] = 9;
	pMatrix[19] = 7;
	pMatrix[20] = 2;
	pMatrix[21] = -1;
	pMatrix[22] = -4;
	pMatrix[23] = 5;
	pMatrix[24] = -7;
}

//  Function  for  memory  allocation  and  definition  of  the  objects elements  
void ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult, double*& x1, double*& x2, int& Size) {
	// Setting the size of the matrix and the vector
	do {
		printf("\nEnter size of the matrix and the vector: ");
		cin >> Size;
		printf("\nChosen size = %d \n", Size);
		if (Size <= 0)
			printf("\nSize of objects must be greater than 0!\n");
	} while (Size <= 0); // Memory allocation
	pMatrix = new double[Size * Size];
	pVector = new double[Size];
	pResult = new double[Size];
	x1 = new double[Size];
	x2 = new double[Size];
	// Initialization of the matrix and the vector elements
	DataInitialization(pMatrix, pVector, Size);
	//RandomDataInitialization(pMatrix, pVector, Size);
	//PRIMER(pMatrix, pVector, Size);

}
// Function for formatted matrix output
void PrintMatrix(double* pMatrix, int RowCount, int ColCount) {
	int i, j; // Loop variables
	for (i = 0; i < RowCount; i++) {
		for (j = 0; j < ColCount; j++)
			printf("%7.4f ", pMatrix[i * RowCount + j]);
		printf("\n");
	}
}
// Function for formatted vector output 
void PrintVector(double* pVector, int Size) {
	int i;
	for (i = 0; i < Size; i++)
		printf("%7.4f ", pVector[i]);
}
// Finding the pivot row
int FindPivotRow(double* pMatrix, int Size, int Iter) {
	int PivotRow = -1; // The index of the pivot row
	int MaxValue = 0; // The value of the pivot element
	int i;             // Loop variable
	// Choose the row, that stores the maximum element
	for (i = 0; i < Size; i++) {
		if ((pSerialPivotIter[i] == -1) &&
			(fabs(pMatrix[i * Size + Iter]) > MaxValue)) {
			PivotRow = i;
			MaxValue = fabs(pMatrix[i * Size + Iter]);
		}
	}
	return PivotRow;
}

// Column elimination 
void SerialColumnElimination(double* pMatrix, double* pVector,
	int Pivot, int Iter, int Size) {
	double PivotValue, PivotFactor;
	PivotValue = pMatrix[Pivot * Size + Iter];
	for (int i = 0; i < Size; i++) {
		if (pSerialPivotIter[i] == -1) {
			PivotFactor = pMatrix[i * Size + Iter] / PivotValue; for (int j = Iter; j < Size; j++) { pMatrix[i * Size + j] -= PivotFactor * pMatrix[Pivot * Size + j]; }pVector[i] -= PivotFactor * pVector[Pivot];
		}
	}
}

// Gaussian elimination
void SerialGaussianElimination(double* pMatrix, double* pVector, int Size) {
	int Iter;       // The number of the iteration of the Gaussian
					// elimination
	int PivotRow;   // The number of the current pivot row
	for (Iter = 0; Iter < Size; Iter++) {
		// Finding the pivot row
		PivotRow = FindPivotRow(pMatrix, Size, Iter);
		pSerialPivotPos[Iter] = PivotRow;
		pSerialPivotIter[PivotRow] = Iter;
		SerialColumnElimination(pMatrix, pVector, PivotRow, Iter, Size);
	}
}

// Back substution 
void SerialBackSubstitution(double* pMatrix, double* pVector,
	double* pResult, int Size) {
	int RowIndex, Row;
	for (int i = Size - 1; i >= 0; i--) {
		RowIndex = pSerialPivotPos[i];
		pResult[i] = pVector[RowIndex] / pMatrix[Size * RowIndex + i];
		for (int j = 0; j < i; j++) {
			Row = pSerialPivotPos[j];
			pVector[Row] -= pMatrix[Row * Size + i] * pResult[i];
			pMatrix[Row * Size + i] = 0;
		}
	}
}
// Function for the execution of Gauss algorithm
void SerialResultCalculation(double* pMatrix, double* pVector,
	double* pResult, int Size) {
	// Memory allocation
	pSerialPivotPos = new int[Size];
	pSerialPivotIter = new int[Size];
	for (int i = 0; i < Size; i++) {
		pSerialPivotIter[i] = -1;
	}
	// Gaussian elimination
	SerialGaussianElimination(pMatrix, pVector, Size);
	// Back substitution
	SerialBackSubstitution(pMatrix, pVector, pResult, Size);
	// Memory deallocation
	delete[] pSerialPivotPos;
	delete[] pSerialPivotIter;
}
// Function for computational process termination
void ProcessTermination(double* pMatrix, double* pVector, double* pResult, double* x1, double* x2) {
	delete[] pMatrix;
	delete[] pVector;
	delete[] pResult;
	delete[] x1;
	delete[] x2;
}
// Function for testing the result 
void TestResult(double* pMatrix, double* pVector,
	double* pResult, int Size) {
	/* Buffer for storing the vector, that is a result of multiplication
	of the linear system matrix by the vector of unknowns */
	double* pRightPartVector;
	// Flag, that shows wheather the right parts
	//vectors are identical or not
	int equal = 0;
	double Accuracy = 1.e-6; // Comparison accuracy 
	pRightPartVector = new double[Size];
	for (int i = 0; i < Size; i++) {
		pRightPartVector[i] = 0;
		for (int j = 0; j < Size; j++) {
			pRightPartVector[i] +=
				pMatrix[i * Size + j] * pResult[j];
		}
	}
	for (int i = 0; i < Size; i++) {
		if (fabs(pRightPartVector[i] - pVector[i]) > Accuracy)
			equal = 1;
	}
	if (equal == 1)
		printf("The result of the parallel Gauss algorithm is NOT correct."
			"Check your code.");
	else
		printf("The result of the parallel Gauss algorithm is correct.")
		; delete[] pRightPartVector;
}

void PrallelZeindelMatrix(double* pMatrix, double* pVector,
	double* pResult, int Size) {
	for (int i = 0; i < Size; i++) {
		for (int j = 0; j < Size; j++) {
			if (i != j)
				pMatrix[i * Size + j] = 1.0 * pMatrix[i * Size + j] / (-pMatrix[i * Size + i]);

		}
		//pVector[i] = pVector[i] / pMatrix[i * Size + i];
	}
	for (int i = 0; i < Size; i++) {
		pVector[i] = 1.0 * pVector[i] / pMatrix[i * Size + i];
	}
	for (int i = 0; i < Size; i++) {
		pMatrix[i * Size + i] = 0.0;
	}
}
double* RaznostVector(double* A, double* B, int Size) {
	for (int i = 0; i < Size; i++) {
		A[i] = A[i] - B[i];
	}
	return A;
}
long double NormaVector(double* A, int Size) {
	long double b = 0.0;
	for (int i = 0; i < Size; i++) {
		b = b + 1.0 * A[i] * A[i];
	}
	return sqrt(b);
}
void Zeidel(double* pMatrix, double* pVector, double* pResult, double* x1, double* x2, int Size) {
	#pragma omp parallel for
	for (int i = 0; i < Size; i++) {
		x1[i] = 0.0;
	}
	while (true) {
		#pragma omp parallel for
		for (int i = 0; i < Size; i++) {
			x2[i] = 0;
		}
		#pragma omp parallel for
		for (int i = 0; i < Size; i++) {
			for (int j = 0; j < Size; j++) {
				x2[i] = x2[i] + pMatrix[i * Size + j] * x1[1];
			}
			x2[i] = x2[i] + pVector[i];
		}
		if (NormaVector(RaznostVector(x1, x2, Size), Size) < epsilon)
			break;
		#pragma omp parallel for
		for (int i = 0; i < Size; i++) {
			x1[i] = x2[i];
		}
	}
	for (int i = 0; i < Size; i++) {
		pResult[i] = x2[i];
	}
}
int main() {
	double* pMatrix;        // The matrix of the linear system
	double* pVector;        // The right parts of the linear system
	double* pResult;        // The result vector
	int     Size;           // The size of the matrix and the vectors
	double start, finish, duration;
	double* x1;
	double* x2;
	// Data initialization
	ProcessInitialization(pMatrix, pVector, pResult, x1, x2, Size);
	start = omp_get_wtime();
	// The matrix and the vectoroutput
		//printf ("Initial Matrix \n");
		//PrintMatrix(pMatrix, Size, Size);
		//printf("Initial Vector \n");
		//PrintVector(pVector, Size);
	//ParallelResultCalculation(pMatrix, pVector, pResult, Size);
	PrallelZeindelMatrix(pMatrix, pVector, pResult, Size);
	Zeidel(pMatrix, pVector, pResult, x1, x2, Size);
	finish = omp_get_wtime();
	duration = finish - start;
	// Testing the result
	//TestResult(pMatrix, pVector, pResult, Size);
	// Printing the time spent by parallel Gauss algorithm
	printf("\n Time of execution: %f\n", duration);
	// The matrix and the vector output
		//printf ("Recalculated Matrix \n");
		//PrintMatrix(pMatrix, Size, Size);
		//printf("Recalculated Vector \n");
		//PrintVector(pVector, Size);
		//printf("\nResult Vector \n");
		//PrintVector(pResult, Size);
		//printf("Parallel Gauss algorithm for solving linear systems\n"); 
	// Program termination
	ProcessTermination(pMatrix, pVector, pResult, x1, x2);
	return 0;
}
