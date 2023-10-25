#include <iostream>
using namespace std;
double** CreateEmptyMatrix(int N){
    double** Matrix = new double * [N];
    for (int i = 0; i <N ; ++i) {
        Matrix[i] = new double[N];
    }
    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            Matrix[i][j] = 0;
        }
    }
    return Matrix;
}
double** CinMatrix(int N){
    double** Matrix = new double * [N];
    for (int i = 0; i <N ; ++i) {
        Matrix[i] = new double[N];
    }
    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
           cin>>Matrix[i][j];
            //Matrix[i][j] = rand()%10;
        }
    }
    return Matrix;
}
void DisplayMatrix(double** Matrix,int N){
    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            cout<<Matrix[i][j]<<"\t";
        }
        cout<<endl;
    }
}
void DeleteMatrix(double** Matrix, int N){

    for (int i = 0; i <N ; ++i) {
        delete[] Matrix[i];
    }
    delete[] Matrix;

}
void MatrixAdd(double ** MatrixA, double** MatrixB, double** MatrixC, int N){

    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            MatrixC[i][j]=MatrixA[i][j]+MatrixB[i][j];
        }
    }
}
void MatrixSubtract(double** MatrixA, double** MatrixB, double** MatrixC, int N){

    for (int i = 0; i <N ; ++i) {
        for (int j = 0; j <N ; ++j) {
            MatrixC[i][j]=MatrixA[i][j]-MatrixB[i][j];
        }
    }
}
void StrassenMultiple(double** MatrixA, double** MatrixB, double** MatrixC, int N){

    double** A11 = CreateEmptyMatrix(N); double** A12 = CreateEmptyMatrix(N); double** A21 = CreateEmptyMatrix(N); double** A22 = CreateEmptyMatrix(N);
    double** B11 = CreateEmptyMatrix(N);double** B12 = CreateEmptyMatrix(N);double** B21 = CreateEmptyMatrix(N);double** B22 = CreateEmptyMatrix(N);
    double** C11 = CreateEmptyMatrix(N);double** C12 = CreateEmptyMatrix(N);double** C21 = CreateEmptyMatrix(N);double** C22 = CreateEmptyMatrix(N);
    double** AA = CreateEmptyMatrix(N); double** BB = CreateEmptyMatrix(N);
    double** M1 = CreateEmptyMatrix(N);
    double** M2 = CreateEmptyMatrix(N);
    double** M3 = CreateEmptyMatrix(N);
    double** M4 = CreateEmptyMatrix(N);
    double** M5 = CreateEmptyMatrix(N);
    double** M6 = CreateEmptyMatrix(N);
    double** M7 = CreateEmptyMatrix(N);

    if(N<3){
        for (int i = 0; i <N ; ++i) {
            for (int j = 0; j <N ; ++j) {
                MatrixC[i][j]=0;
                for (int k = 0; k <N ; ++k) {
                    MatrixC[i][j]=MatrixC[i][j] + MatrixA[i][k] * MatrixB[k][j];
                }
            }
        }
    }
    else{
        for (int i = 0; i <N/2 ; ++i) {
            for (int j = 0; j <N/2 ; ++j) {
                A11[i][j]=MatrixA[i][j];
                A12[i][j]=MatrixA[i][j + N/2];
                A21[i][j]=MatrixA[i + N/2][j];
                A22[i][j]=MatrixA[i + N/2][j + N/2];

                B11[i][j]=MatrixB[i][j];
                B12[i][j]=MatrixB[i][j + N/2];
                B21[i][j]=MatrixB[i + N/2][j];
                B22[i][j]=MatrixB[i + N/2][j + N/2];
            }
        }


        //обчислення M1 = (A11 + A22) × (B11 + B22)
        MatrixAdd(A11, A22, AA, N/2);
        MatrixAdd( B11, B22, BB, N/2);
        StrassenMultiple( AA, BB, M1, N/2);

        //обчислення M2 = (A21 + A22) × B11
        MatrixAdd( A21, A22, AA, N/2);
        StrassenMultiple( AA, B11, M2, N/2);

        //обчислення M3 = A11 × (B12 - B22)
        MatrixSubtract( B12, B22, BB, N/2);
        StrassenMultiple( A11, BB, M3, N/2);

        //обчислення M4 = A22 × (B21 - B11)
        MatrixSubtract( B21, B11, BB, N/2);
        StrassenMultiple( A22, BB, M4, N/2);

        //обчислення M5 = (A11 + A12) × B22
        MatrixAdd( A11, A12, AA, N/2);
        StrassenMultiple( AA, B22, M5, N/2);

        //обчислення M6 = (A21 - A11) × (B11 + B12)
        MatrixSubtract( A21, A11, AA, N/2);
        MatrixAdd( B11, B12, BB, N/2);
        StrassenMultiple( AA, BB, M6, N/2);

        //обчислення M7 = (A12 - A22) × (B21 + B22)
        MatrixSubtract( A12, A22, AA, N/2);
        MatrixAdd( B21, B22, BB, N/2);
        StrassenMultiple( AA, BB, M7, N/2);

        //обчислення C11 = M1 + M4 - M5 + M7
        MatrixAdd( M1, M4, AA, N/2);
        MatrixSubtract( AA , M5, BB, N/2);
        MatrixAdd( BB, M7, C11, N/2);

        //обчислення C12 = M3 + M5
        MatrixAdd( M3, M5, C12, N/2);

        //обчислення C21 = M2 + M4
        MatrixAdd(M2, M4, C21, N/2);

        //обчислення C22 = M1 - M2 + M3 + M6
        MatrixSubtract( M1, M2, AA, N/2);
        MatrixAdd( AA, M3, BB, N/2);
        MatrixAdd( BB, M6, C22, N/2);

        //обчислення MatrixС[N][N]
        for (int i = 0; i < N / 2; i++) {
            for (int j = 0; j < N / 2; j++) {
                MatrixC[i][j] = C11[i][j];
                MatrixC[i][j + N / 2] = C12[i][j];
                MatrixC[i + N / 2][j] = C21[i][j];
                MatrixC[i + N / 2][j + N / 2] = C22[i][j];
            }
        }
        for (int i = 0; i < N; ++i) {
            delete[] A11[i];
            delete[] A12[i];
            delete[] A21[i];
            delete[] A22[i];

            delete[] B11[i];
            delete[] B12[i];
            delete[] B21[i];
            delete[] B22[i];

            delete[] C11[i];
            delete[] C12[i];
            delete[] C21[i];
            delete[] C22[i];

            delete[] AA[i];
            delete[] BB[i];

            delete[] M1[i];
            delete[] M2[i];
            delete[] M3[i];
            delete[] M4[i];
            delete[] M5[i];
            delete[] M6[i];
            delete[] M7[i];
        }
    }
}
int main(){
    int N = 4;
    cout<<"enter matrix size"<<endl;
    cout<<"Matrix A = "<<endl;
    double** MatrixA = CinMatrix(N);
    DisplayMatrix(MatrixA, N);
    cout<<endl;
    cout<<"Matrix B = "<<endl;
    double** MatrixB = CinMatrix(N);
    DisplayMatrix(MatrixB, N);
    double** MatrixC = CreateEmptyMatrix(N);
    cout<<"the result of strassen`s algorithm Matrix C = "<<endl;
    StrassenMultiple(MatrixA, MatrixB, MatrixC, N);
    DisplayMatrix(MatrixC, N);
    DeleteMatrix(MatrixA, N);
    DeleteMatrix(MatrixB, N);
    DeleteMatrix(MatrixC, N);

    return 0;
}