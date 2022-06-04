#include "mex.h"
#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <limits>  
using namespace std;

#define A_IN prhs[0]
#define B_IN prhs[1]
#define C_IN prhs[2]
#define x1_OUT plhs[0]
#define x2_OUT plhs[1]
#define x3_OUT plhs[2]
#define x4_OUT plhs[3]

void myprint(double *A, size_t n) {
    mexPrintf("matrix\n");
    for (size_t i = 0; i < n; i++) {
        char buffer[20];
        size_t k=sprintf (buffer, "%f |", A[i]);
        mexPrintf(buffer, k);
    }
    mexPrintf("\n");
    return;
}

#define myassert(isOK, astr)                                            \
  do {                                                                  \
    if(!(isOk)) {                                                       \
      std::ostringstream fmt;                                           \
      fmt << "In " << __PRETTY_FUNCTION__ << ", "                       \
          << __FILE__ << ":" << __LINE__ << ": " << (astr);             \
      (void) mexErrMsgTxt(fmt.str().c_str());                           \
    }                                                                   \
  } while(false)

void biquadsolve(
        vector<complex<double>>& A, 
        vector<complex<double>>& B, 
        vector<complex<double>>& C, 
        double *outMatrix_x1_real, 
        double *outMatrix_x1_imag, 
        double *outMatrix_x2_real, 
        double *outMatrix_x2_imag, 
        double *outMatrix_x3_real, 
        double *outMatrix_x3_imag,
        double *outMatrix_x4_real, 
        double *outMatrix_x4_imag) 
{
    for (size_t i = 0; i < A.size(); i++) {
            complex<double> a = A[i], b = B[i], c = C[i];
            complex<double> d = b * b - 4.0 * a * c;
            complex<double> x1, x2, x3, x4, tmp1, tmp2;
            double epsilon = 1e-5;// for doubles comparison
            if (abs(a) > epsilon) {
                tmp1 = (-b - sqrt(d)) / (2.0 * a);
                tmp2 = (-b + sqrt(d)) / (2.0 * a);
                x1 = sqrt(tmp1);
                x2 = -x1;
                x3 = sqrt(tmp2);
                x4 = -x3;
            }
            else if (abs(b) > epsilon) {
                x1 = sqrt(-c / b);
                x2 = -x1;
                x3 = numeric_limits<double>::quiet_NaN();
                x4 = numeric_limits<double>::quiet_NaN();
            }
            else {
                x1 = numeric_limits<double>::quiet_NaN();
                x2 = numeric_limits<double>::quiet_NaN();
                x3 = numeric_limits<double>::quiet_NaN();
                x4 = numeric_limits<double>::quiet_NaN();
            }
            outMatrix_x1_real[i] = x1.real();
            outMatrix_x1_imag[i] = x1.imag();
            outMatrix_x2_real[i] = x2.real();
            outMatrix_x2_imag[i] = x2.imag();
            outMatrix_x3_real[i] = x3.real();
            outMatrix_x3_imag[i] = x3.imag();
            outMatrix_x4_real[i] = x4.real();
            outMatrix_x4_imag[i] = x4.imag();
        }
}

void mexFunction(
        int nlhs, 
        mxArray *plhs[],
        int nrhs,
        const mxArray *prhs[])
{
    char buffer [100];
    int n;
    size_t nrows, ncols;                   /* size of matrices */
    
    if(nrhs!=3) { 
        mexErrMsgIdAndTxt("MyToolbox:biquadsolve:nrhs","Three inputs required.");
    }
    
    if(nlhs != 2 && nlhs != 4) {
        mexErrMsgIdAndTxt("MyToolbox:biquadsolve:nlhs","Must have either 2 or 4 output arguments.");
    }

    if( (mxGetM(prhs[0])!= mxGetM(prhs[1])) || (mxGetM(prhs[1])!= mxGetM(prhs[2]))
    ||(mxGetN(prhs[0])!= mxGetN(prhs[1])) || (mxGetN(prhs[1])!= mxGetN(prhs[2])) ) {
        mexErrMsgIdAndTxt("MyToolbox:biquadsolve:nrhs","Input matrices must be the same size.");
    }
     /* get dimensions of the input matrix */
    nrows = mxGetM(A_IN);
    ncols = mxGetN(A_IN);
    double *A1;
    A1 = mxGetPr(prhs[1]);
    

    double *real_data_ptr_A = (double *)mxGetPr(A_IN);
    double *imag_data_ptr_A = (double *)mxGetPi(A_IN);
    double *real_data_ptr_B = (double *)mxGetPr(B_IN);
    double *imag_data_ptr_B = (double *)mxGetPi(B_IN);
    double *real_data_ptr_C = (double *)mxGetPr(C_IN);
    double *imag_data_ptr_C = (double *)mxGetPi(C_IN);
    vector<complex<double>> A, B, C;
    double *outMatrix_x1_real, *outMatrix_x1_imag, *outMatrix_x2_real, *outMatrix_x2_imag, *outMatrix_x3_real, *outMatrix_x3_imag, *outMatrix_x4_real, *outMatrix_x4_imag;
    
    if( !mxIsComplex(prhs[0]) || !mxIsComplex(prhs[1]) || !mxIsComplex(prhs[2]) ) {
        //mexPrintf("not a complex\n");
        for (size_t i = 0; i < nrows*ncols; i++)
            {
                A.push_back(complex<double>(real_data_ptr_A[i], 0));
                B.push_back(complex<double>(real_data_ptr_B[i], 0));
                C.push_back(complex<double>(real_data_ptr_C[i], 0));
            }
        
        //mexErrMsgIdAndTxt( "MyToolbox:quadsolve:nrhs", "Inputs must be complex.\n");
    } else { 
        for (size_t i = 0; i < nrows*ncols; i++)
            {
                A.push_back(complex<double>(real_data_ptr_A[i], imag_data_ptr_A[i]));
                B.push_back(complex<double>(real_data_ptr_B[i], imag_data_ptr_B[i]));
                C.push_back(complex<double>(real_data_ptr_C[i], imag_data_ptr_C[i]));

            }
    }
    plhs[0] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxCOMPLEX);
    plhs[2] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxCOMPLEX);
    plhs[3] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxCOMPLEX);
    outMatrix_x1_real = mxGetPr(x1_OUT);
    outMatrix_x1_imag = mxGetPi(x1_OUT);
    outMatrix_x2_real = mxGetPr(x2_OUT);
    outMatrix_x2_imag = mxGetPi(x2_OUT);
    outMatrix_x3_real = mxGetPr(x3_OUT);
    outMatrix_x3_imag = mxGetPi(x3_OUT);
    outMatrix_x4_real = mxGetPr(x4_OUT);
    outMatrix_x4_imag = mxGetPi(x4_OUT);
   
    biquadsolve(A, B, C, outMatrix_x1_real, outMatrix_x1_imag, outMatrix_x2_real, outMatrix_x2_imag, outMatrix_x3_real, outMatrix_x3_imag, outMatrix_x4_real, outMatrix_x4_imag);

return;
}