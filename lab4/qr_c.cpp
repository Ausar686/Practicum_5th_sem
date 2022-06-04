#include "mex.h"
#include "matrix.h"
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>

using namespace std;


#define A_IN prhs[0]
#define Q_OUT plhs[0]
#define R_OUT plhs[1]
#define eps 1/1000

int getrk(vector<double> A){
    int n = int(sqrt(A.size()));
    int rg = n;
    for (int i = 0; i < n; i++){
        if (!A[i*n + i]){
            int j;
            for (j = i+1; j < n && !A[j*n+i]; j++);
                if (j == n){
                    rg--;
                    continue;
                }
                else
                    for (int k = i; k < n; k++){
                        bool t = A[i*n+k];
                        A[i*n+k] = A[j*n+k];
                        A[j*n+k] = t;
                    }
        }
        for (int j = i+1; j < n; j++){
            if (A[j*n+i]){
                for (int k = i; k < n; k++)
                    A[j*n+k] = A[j*n+k] - A[i*n+k];
            }
        }
    }
    return rg;
}

double norm(vector<double> X)
{
    double s = 0;
    int n = X.size();
    for(int i=0;i<n;i++)
    {
        s = s + X[i] * X[i];
    }
    s = sqrt(s);
    return s;
}
double scl(vector<double> X, vector<double> Z)
{
    double s = 0;
    int n = X.size();
    for(int i=0;i<n;i++)
    {
        s = s+X[i]*Z[i];
    }
    return s;
}
vector<double> col( vector<vector<double>> &A,int i)
{
    vector<double> res;
    int j;
    int n = A.size();
    for(j=0; j<n; j++)
    {
        res.push_back(A[j][i]);
    }
    return res;
}
void qr_c(vector<double>& A, double *Q, double *R)
{
    int n = int(sqrt(A.size())); // квадратная матрица, поэтому можем просто взять корень
    vector<double> c(n,0);
    vector<vector<double>> a,b;
    int i,j,k;
    for(i=0; i<n; i++)
    {
        a.push_back(vector<double>());
        b.push_back(vector<double>()); // кладём пустой вектор
        for(j=0; j<n; j++)
        {
            b[i].push_back(A[i+n*j]);
            a[i].push_back(A[i+n*j]);
        }
    }
    double norma0 = norm(col(b,0));
    for(i=0; i<n; i++){
        b[i][0] = b[i][0]/norma0;
    }
    i = 1;
    while(i<n)
    {
        for(k=0; k<i; k++)
        {
            c[k] = scl(col(b,k),col(b,i))/scl(col(b,k), col(b,k));
            for(j=0; j<n; j++){
                b[j][i] -= c[k]*b[j][k];
            }
        }
        double norma_i = norm(col(b,i));
        for(j=0; j<n; j++){
            b[j][i] /= norma_i;
        }
        i++;
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            Q[i+n*j] = -b[i][j];
        }
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            double s = 0;
            vector<double> col1 = col(b,i);
            vector<double> col2 = col(a,j);
            for(k=0; k<n;k++)
            {
               s += col1[k]*col2[k];
            }
            R[i+n*j] = -s;
        }
    }
    int flag = 1;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(abs(A[i+n*j]) > eps)
                flag = 0;
        }
    }
    if ((flag==0)&&(getrk(A) < n))
    {
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                Q[i+n*j] = 0;
                R[i+n*j] = A[i+n*j];
            }
            Q[i+n*i] = 1;
        }
    } else if (flag==1)
    {
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                R[i+n*j] = 0;
                Q[i+n*j] = 0;
            }
            Q[i+n*i] = 1;
        }

    }
    return ;
}
void mexFunction(int nlhs, mxArray *plhs[],int nrhs,  const mxArray *prhs[])
{
    size_t nrows, ncols;
    // Проверка корректности
    if(nrhs!=1)
        mexErrMsgIdAndTxt("MyToolbox:qr_c:nrhs","Неверное число аргументов.");
    if(mxGetM(prhs[0])!= mxGetN(prhs[0]))
        mexErrMsgIdAndTxt("MyToolbox:qr_c:nrhs","Матрица не квадратная!");
    if(!mxIsNumeric(prhs[0])) 
        mexErrMsgIdAndTxt( "MyToolbox:qr_c:plhs", "Матрица заполнена некорректно!");

    nrows = mxGetM(A_IN);
    ncols = mxGetN(A_IN);
    vector<double> B;

    double *A = (double *)mxGetPr(A_IN);
    double *Q, *R;
    for (size_t i = 0; i < nrows*ncols; i++)
        {
            B.push_back(A[i]);
        }
    plhs[0] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxREAL);
    Q = mxGetPr(Q_OUT);
    R = mxGetPr(R_OUT);
    qr_c(B, Q, R);
return ;
}