#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;
using namespace arma;

void print_matrix(double **M, int rows, int cols){
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
}

void delete_matrix(double **M, int rows){
    for (int i = 0; i < rows; i++){
        delete[] M[i];
    }
    delete[] M;

}

void tridiagonal_same(double **M, int dimension, double diagonal,
                        double under_diag, double over_diag){
    for(int i = 0; i < dimension-1; i++){
        M[i][i] = diagonal;
        M[i+1][i] = under_diag;
        M[i][i+1] = over_diag;
    }
    M[dimension-1][dimension - 1] = diagonal;
}

mat tridiagonal_same_armadillo(int size, double diagonal,
                        double lower_diagonal, double upper_diagonal){
    mat A(size, size, fill::zeros);
    A.eye();
    A = diagonal*A;
    for(int i = 1; i < size; i++){
        A(i-1,i) = upper_diagonal;
        A(i, i-1) = lower_diagonal;
    }
    return A;
}

//The source function used in this project
double myFunction(double x){
    return 100*exp(-10*x);
}
//The exact solution of the second order differential equation
double exact_double_derivative(double x){
    return 1- (1-exp(-10))*x-exp(-10*x);
}

//Simple algorithm for solving a tridiagonal matrix with -1's on
//the lower and upper diagonal and two-point boundary condition where
//the endpoints are 0 (Dirichlet boundary condition).
void simple_triDiag_solver(double *a, double *g, double *v, int length){
    //Forward substitution
    for(int i = 2; i < length; i++){
        a[i] = a[i] - 1/a[i-1];
        g[i] = g[i] + g[i-1]/a[i-1];
    }
    //Backward substitution
    v[length-2] = g[length-2]/a[length-2];
    for(int i = length-3; i > 0; i--){
        v[i] = (g[i]+v[i+1])/(a[i]);
    }
}

int main()
{
    //Variables for performancetesting
    clock_t start, finish;
    //Initializing needed parameters
    int length_number_of_points = 7;
    int number_of_points_array[length_number_of_points] =
    {10, 100, 1000, 10000, 100000, 1000000, 10000000};
//    string number[5] = {"10", "100", "1000", "10000", "100000"};
    double errors_for_different_step_lengths[length_number_of_points];
    double step_lengths[length_number_of_points];
    cout << "Time my algo\t" << "Time LU armadillo" << endl;
    for(int j = 0; j < length_number_of_points; j++){
        int number_of_points = number_of_points_array[j];
        double step_length = 1./(number_of_points-1);
        double diagonal = 2;
        double lower_diagonal = -1;
        double upper_diagonal = -1;

        //My algorithm
        double *diagonal_elements;
        diagonal_elements = new double[number_of_points];
        for(int i = 0; i < number_of_points; i++){
            diagonal_elements[i] = diagonal;
        }
        double *source_function,*v,*u_exact,*u_estimated,*error;
        source_function = new double[number_of_points];
        v = new double[number_of_points];
        u_exact = new double[number_of_points];
        u_estimated = new double[number_of_points];
        error = new double[number_of_points];
        //    Initializing source functions for my algorithm and armadillo
        for(int i = 0; i < number_of_points; i++){
            source_function[i] = myFunction(i*step_length);
            u_exact[i] = exact_double_derivative(i*step_length);
            v[i] = 0;
        }
        vec source_armadillo(number_of_points-2, fill::zeros);
        for(int i = 0; i < number_of_points-2; i++){
            source_armadillo(i) = myFunction((i+1)*step_length);
        }

        //Solving the problem with my algorithm
        start = clock();
        simple_triDiag_solver(diagonal_elements, source_function,
                              v, number_of_points);
        finish = clock();
        cout << setprecision(7) << (finish-start)/double(CLOCKS_PER_SEC) << "\t\t";
        for(int i = 0; i < number_of_points; i++){
            u_estimated[i] = step_length*step_length*v[i];
            error[i] = log10(abs((u_estimated[i]-u_exact[i])/u_exact[i]));
        }
        error[number_of_points-1] = 0;
        error[0] = 0;
        double error_max = -100000;
        for(int i = 1; i < number_of_points-1; i++){
            if(error[i]>error_max){
                error_max = error[i];
            }
        }
        errors_for_different_step_lengths[j] = error_max;
        step_lengths[j] = step_length;

        if(number_of_points <= 10000){
            mat A = tridiagonal_same_armadillo(number_of_points-2, diagonal,
                                               lower_diagonal, upper_diagonal);
            vec x = solve(A, source_armadillo);
            vec v_estimated_armadillo(number_of_points, fill::zeros);
            for(int i = 1; i < number_of_points-1; i++){
                v_estimated_armadillo(i) = x(i-1);
            }
            //Solving the problem with LU decomposition
            mat L, U;
            lu(L, U, A);
            start = clock();
            vec tmp = solve(L, source_armadillo);
            vec v_arma_lu = solve(U,tmp);
            finish = clock();
            cout << right << setprecision(7) << (finish-start)/double(CLOCKS_PER_SEC) << endl;
            vec v_lu_estimate(number_of_points, fill::zeros);
            for(int i = 1; i < number_of_points-1; i++){
                v_lu_estimate(i) = v_arma_lu(i-1);
            }
        }else{
            cout << right << "Not enough memory for matrix" << endl;
        }
        string s = "exact_and_approximated_n"+to_string(number_of_points_array[j])+".txt";
        //Writing results to file
        ofstream myfile;
        myfile.open(s.c_str());
        myfile << number_of_points << "\t" << errors_for_different_step_lengths[j] << endl;
        for(int i = 0; i < number_of_points; i++){
            myfile << right << setw(26) << setprecision(5) << u_exact[i] <<
                      "\t"<< u_estimated[i] << "\t" << error[i] << endl;
        }
        myfile.close();
    }
    ofstream myfile;
    myfile.open("errors_for_different_steplengths.txt");
    myfile << length_number_of_points << endl;
    for(int i = 0; i < length_number_of_points; i++){
        myfile << right << setw(26) << setprecision(5) << step_lengths[i] << "\t" <<
                  errors_for_different_step_lengths[i] << endl;
    }
    myfile.close();
}





