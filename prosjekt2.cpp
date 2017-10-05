#include <iostream>
#include <cmath>   
#include <armadillo>  
#include <assert.h>

using namespace arma;

//Finding maximum element of off-diagonals
double maxoffdiag( mat &A, int *k, int *l, int n )
{
    double max = 0.0;

    for(int i = 0; i < n; i++){
        for( int j = i + 1; j < n; j++){
            if ( fabs(A(i, j)) > max ){
                max = fabs(A(i, j));
                *l = i;
                *k = j;
            }
        }
    }
    return max;
} 

void rotate ( mat &A, mat &R, int k, int l, int n)  //TODO why no pointers at l and k?
{
    double s, c;
    if ( A(k, l) != 0.0 ){
        double t, tau;
        tau = (A(l,l) - A(k, k))/(2*A(k, l));
        if (tau > 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        } else {
            t = -1.0/( -tau + sqrt(1.0 + tau*tau));
        }                                  

        c= 1/sqrt(1+t*t);
        s = c*t; 

    } else {
        c = 1.0;
        s = 0.0;
    }

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il; //Because we want to be able to change elements and still be able to use the previous one in the next step. 
    a_kk = A(k,k);
    a_ll = A(l,l);
    //changing the matric elements with indices k and l, maximum off diagonal elements?
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0; //hard coding of the zeros, yay
    A(l,k) = 0.0; //symetric matrix

    //change the remaining elements
    for ( int i = 0; i < n; i++ ){
        if( i != k && i != l ){
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l); 
        }
        //compute new eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return;
}
    

void jacobi_method (mat &A, mat &R, int n)
{
    // making identity matrix
    R.eye();

    int k, l;
    double epsilon = 1.0e-8;
    double max_number_iterations = (double) n * (double) n * (double) n;
    int iterations = 0;
    double max_offdiag = maxoffdiag(A, &k, &l, n);

    while (fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations){
        max_offdiag = maxoffdiag(A, &k, &l, n);
        rotate(A, R, k, l, n);
        iterations++;
    }
    std::cout << "Number of iterations; " <<iterations << "\n";

    
    return;
}   

vec sort_eigen(mat A, int n){ 
    vec eigenvalues = zeros<vec>(n); 
    eigenvalues = sort(A.diag());  

    return eigenvalues;
}

//Test functions
//----------------------------------------------------------------------------------

void testMax(){
    // takes a known matrix with known max value on off-diagonal. Checks that maxoffdiag gives correct answer. 
    int n = 5;
    int k,l;
    mat test = zeros(n, n);

    for( int i = 0; i < n; i++){
        for( int j = 0; j < n; j++){
            test(i, j) = i + j;
        }
    } 

    double maxTest = maxoffdiag(test, &k, &l, n);

    assert (maxTest == 7);
    std::cout <<"maxtest passed" << endl;
    
    return;
}

void testEigenvalues(){
    int n = 2;
    double epsilon = 1E-8;
    mat testR = zeros(n,n);
    mat testA = zeros(n,n);

    testA(0,0) = 1.;
    testA(0,1) = 4.;
    testA(1,0) = 4.;
    testA(1,1) = 7.;

    jacobi_method(testA, testR, n); 
    vec eigenvalues = sort_eigen(testA, n);    
    vec analytical_eigen = zeros<vec>(n);

    analytical_eigen(0) = -1;
    analytical_eigen(1) = 9;

    for( int i = 0; i < n; i++){
        assert (std::abs(eigenvalues(i) - analytical_eigen(i)) <= epsilon); 
    }
    std::cout <<"eigenvalue test passed" << endl;


//Main
//-------------------------------------------------------------------------------

}

int main( int argc, char *argv[]){
    
    int N = atoi(argv[1]);

    double rho_min = 0.;
    double rho_max = 5.;
    double h = (rho_max - rho_min)/((double) N+1); 

    //std::cout << h << endl;



    // vec rho = zeros<vec>(N);   
    mat A = zeros(N, N);
    vec d = zeros<vec>(N);
    mat R = zeros(N, N);

    for( int i = 0; i < N; i++){ 
        double rho = rho_min + (i+1)*h;
        d(i) = 2/(h*h) + rho*rho;
    }


    A.diag(-1) += -1/(h*h); 
    A.diag(1) += -1/(h*h);
    A.diag(0) = d;

    testMax(); 
    testEigenvalues();
    jacobi_method(A, R, N);

    
    //std::cout << eigenvalues << endl; 
}

                 
