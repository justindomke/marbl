#include <stdio.h>
#include <map>
#include <Eigen/Eigen>
//#include <Eigen/Array>
#include <ctime>
#include <iostream>
#include <cfloat>

using namespace std;
using namespace Eigen;

typedef Map<MatrixXd> MatrixMd;
typedef Map<MatrixXi> MatrixMi;


// intended for eigen matrix types
template <typename T>
void norm_cols(T& A){
//     cout << "A:\n" << A << endl;
//     for(int i=0; i<A.cols(); i++){
//         double mysum = 0;
//         for(int j=0; j<A.rows(); j++){
//             cout << "i: " << i << " j: " << j << " A[j,i]: " << A[j,i] << endl;
//             mysum += A[j,i];
//         }
//         cout << "mysum: " << mysum << endl;
//         for(int j=0; j<A.rows(); j++)
//             A[j,i] /= mysum;
//     }
    for(int i=0; i<A.cols(); i++)
        A.col(i) = A.col(i) / A.col(i).sum();
}

// actually, this should work for any matrix type!
template <typename T>
T imresize_bigger(const T& A){
    int M  = A.rows();
    int N  = A.cols();
    int M2 = M*2;
    int N2 = N*2;
    
        // want f(M2-1)=M-1-r
        //      f(0)   = r
        //      f(x)   = a*x+b
        //      b      = r;
        //      a*(M2-1)+r = M-1-r
        //      a          = (M-1-2*r)/(M2-1)
    
    T C = T(M2,N2);
    for(int y2=0; y2<M2; y2++){
        for(int x2=0; x2<N2; x2++){
            double r = -.25;
            double a = (M-1-2*r)/(M2-1);
            double b = r;
            double y = a*y2+b;
            a = (N-1-2*r)/(N2-1);
            b = r;
            double x = a*x2+b;
            if(y<0) y=0;
            if(x<0) x=0;
            if(y>M-1) y=M-1;
            if(x>N-1) x=N-1;
            int y_lo = floor(y);
            int x_lo = floor(x);
            int y_hi = ceil(y);
            int x_hi = ceil(x);
            double yfrac = 1-(y-y_lo);
            double xfrac = 1-(x-x_lo);
            
            C(y2, x2) =    yfrac *   xfrac *A(y_lo,x_lo) + 
                        (1-yfrac)*   xfrac *A(y_hi,x_lo)  + 
                           yfrac *(1-xfrac)*A(y_lo,x_hi)  + 
                        (1-yfrac)*(1-xfrac)*A(y_hi,x_hi);
        }
    }
    
    return C;
}

template <typename T>
T imresize_smaller(T A){
    int M = A.rows();
    int N = A.cols();
    int M2 = M/2;
    int N2 = N/2;
    
    //cout << "M2: " << M2 << "  N2: " << N2 << endl;
    
    T C = T(M2,N2);
    
    for(int y2=0; y2<M2; y2++){
        for(int x2=0; x2<N2; x2++){
            int y_lo = y2*2;
            int x_lo = x2*2;
            int y_hi = y_lo+1;
            int x_hi = x_lo+1;
            C(y2, x2) = .25*A(y_lo,x_lo) + 
                        .25*A(y_hi,x_lo) + 
                        .25*A(y_lo,x_hi) + 
                        .25*A(y_hi,x_hi);
        }
    }
    return C;
}

template <typename T>
T rescale(T A, int diff){
    if(diff==0)
        return A;
    if(diff<0) // make bigger
        return rescale(imresize_bigger(A),diff+1);
    else
        return rescale(imresize_smaller(A),diff-1);
}

// templated function to print an arbitrary list of crap
template <typename T>
void printStuff(T first, T last){
     for(; first != last; ++first)
         cout << *first << endl;
} //use like: printStuff(ingredients.begin(), ingredients.end()

// // this is done in a somewhat stupid way returning a vector just because
// // I don't feel like remembering how to use templates properly
// template <typename T>
// T log_sum_exp(T A){
//     T rez = T(1,1);
//     T damin = T(1,1);
//     damin(0) = A.minCoeff();
//     rez(0) = 0;
//     for(int i=0; i<A.rows()*A.cols(); i++)
//         rez(0) += exp(A(i)-damin(0));
//     rez(0) = log(rez(0))+damin(0);
//     return rez;
// }

// previously, this looked for the minimum number, which seems foolish...
inline double log_sum_exp(const MatrixXd& A){
    int nvals = A.rows()*A.cols();
    
    double rez   = 0;
    //double damin = A.maxCoeff();
    double damin = DBL_MIN;
    for(int i=0; i<nvals; i++)
        damin = max(damin,A(i));
    
    for(int i=0; i<nvals; i++)
        rez += exp(A(i)-damin);
    rez = log(rez)+damin;
    return rez;    
}

inline double log_sum_exp(const MatrixXd& A, int nvals){
    double rez   = 0;
    //double damin = A.maxCoeff();
    double damin = DBL_MIN;
    for(int i=0; i<nvals; i++)
        damin = max(damin,A(i));
    
    for(int i=0; i<nvals; i++)
        rez += exp(A(i)-damin);
    rez = log(rez)+damin;
    return rez;    
}


// void print_tapestats(){
//     int counts[11];
//     tapestats(1, counts);
//     cout << "#ind:  " << counts[0] << endl;
//     cout << "#dep:  " << counts[1] << endl;
//     cout << "#live: " << counts[2] << endl;
//     cout << "#vstk: " << counts[3] << endl;
//     cout << "#bsiz: " << counts[4] << endl;
//     cout << "#ops:  " << counts[5] << endl;
//     cout << "#6:    " << counts[6] << endl;
//     cout << "#7:    " << counts[7] << endl;
//     cout << "#8:    " << counts[8] << endl;
//     cout << "#9:    " << counts[9] << endl;
//     cout << "#10:   " << counts[10] << endl;    
// }

// template <typename T>
// vector<MatrixXd> niceGrad(T first, T last){
//     // first, count the total number of elements
//     int numels = 0;
//     for(T it=first; it != last; ++it)
//         numels += it->rows()*it->cols();
//     //cout << "numels: " << numels << endl;
//     // declare a big vector to stick everything in
//     double x[numels];
//     double g[numels];
//     int el=0;
//     for(T it=first; it != last; ++it)
//         for(int i=0; i<it->rows()*it->cols(); i++)
//             x[el++] = (*it)(i);
//     // call gradient function
//     //double f;
//     //function(1,1,numels,x,&f);  
//     clock_t start = clock();
//     gradient(1  ,numels,x,g);
//     cout << "time3: " << (double(clock())-double(start))/CLOCKS_PER_SEC << endl;
//     //for(int i=0; i<numels; i++)
//     //    cout << x[i] << " " << g[i] << endl;
//     // create output matrices holding gradient values
//     vector<MatrixXd> vec;
//     el = 0;
//     for(T it=first; it != last; ++it){
//         // new matrix of same size
//         MatrixXd M(it->rows(),it->cols());
//         for(int i=0; i<it->rows()*it->cols(); i++)
//             M(i) = g[el++];
//         vec.push_back(M);
//     }
//     return vec;
// }
// 
// template <typename T>
// MatrixXa adMat(T& x0){
//     //MatrixXa x = x0.cast<adouble>(); // not sure why doesn't work
//     MatrixXa x(x0.rows(),x0.cols());
//     for(int i=0; i<x.rows()*x.cols(); i++){
//         x[i] <<= x0[i];
//         //x[i] = x0[i];
//         //x[i].declareIndependent();
//     }
//     return x;
// }
