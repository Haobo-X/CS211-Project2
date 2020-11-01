#include "../include/for_you_to_do.h"
#include <math.h>

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 126;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    int i, j, k, max_index, tmp2;
    double max, tmp1;
    // used for swapping rows
    register double * tmp_row = (double *)malloc(sizeof(double) * n);
    
    for (i = 0; i < n; i++)
    {
        max_index = i;
        max = fabs(A[i * n + i]);
        for (j = i + 1; j < n; j++)
        {
            tmp1 = fabs(A[j * n + i]);
            if (max < tmp1)
            {
                max_index = j;
                max = tmp1;
            }
        }
        
        // if the matrix is singular
        if (max == 0)
        {
            return -1;
        }
            
        if (max_index != i)
        {
            tmp2 = ipiv[i];
            ipiv[i] = ipiv[max_index];
            ipiv[max_index] = tmp2;
            // swap rows of A
            memcpy(tmp_row, A + i * n, n * sizeof(double));
            memcpy(A + i * n, A + max_index * n, n * sizeof(double));
            memcpy(A + max_index * n, tmp_row, n * sizeof(double));
        }

        for (j = i + 1; j < n; j++)
        {
            A[j * n + i] = A[j * n + i] / A[i * n + i];
            for (k = i + 1; k < n; k++)
            {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];
            }
        }
    }
    free(tmp_row);
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
    int i, j;
    register double sum = 0;
    register double * tmp_B = (double *)malloc(sizeof(double) * n);
    
    if (UPLO == 'L')
    {
        for (i = 0; i < n; i++)
        {
            tmp_B[i] = B[ipiv[i]];
        }
        
        for (i = 0; i < n; i++)
        {
            sum = tmp_B[i];
            for (j = 0; j < i; j++)
            {
                sum -= B[j] * A[i * n + j];
            }
            B[i] = sum;
        }
    }
    else if (UPLO == 'U')
    {
        for (i = n - 1; i >= 0; i--)
        {
            sum = 0;
            for (j = i + 1; j < n; j++)
            {
                sum += B[j] * A[i * n + j];
            }
            B[i] = (B[i] - sum) / A[i * n + i];
        }
    }
    free(tmp_B);
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
/*
    int i1, j1, k1, i2, j2, k2;
    for (i1 = 0; i1 < n; i1 += b)
    {
        for (j1 = 0; j1 < j; j1 += b)
        {
            for (k1 = 0; k1 < i; k1 += b)
            {
                for (i2 = i1; i2 < (i1 + b > n? n : (i1 + b)); i2 += 3)
                {
                    for (j2 = j1; j2 < (j1 + b > j? j : (j1 + b)); j2 += 3)
                    {
                        register double C_0_0 = C[i2 * k + j2];
                        register double C_0_1 = C[i2 * k + (j2 + 1)];
                        register double C_0_2 = C[i2 * k + (j2 + 2)];
                        register double C_1_0 = C[(i2 + 1) * k + j2];
                        register double C_1_1 = C[(i2 + 1) * k + (j2 + 1)];
                        register double C_1_2 = C[(i2 + 1) * k + (j2 + 2)];
                        register double C_2_0 = C[(i2 + 2) * k + j2];                
                        register double C_2_1 = C[(i2 + 2) * k + (j2 + 1)];
                        register double C_2_2 = C[(i2 + 2) * k + (j2 + 2)];

                        for (k2 = k1; k2 < (k1 + b > i? i : (k1 + b)); k2++)
                        {
                            register double A_0 = A[i2 * k + k2];
                            register double A_1 = A[(i2 + 1) * k + k2];
                            register double A_2 = A[(i2 + 2) * k + k2];
                            register double B_0 = B[k2 * k + j2];
                            register double B_1 = B[k2 * k + (j2 + 1)];
                            register double B_2 = B[k2 * k + (j2 + 2)];
                            
                            C_0_0 += A_0 * B_0;
                            C_0_1 += A_0 * B_1;
                            C_0_2 += A_0 * B_2;
                            C_1_0 += A_1 * B_0;
                            C_1_1 += A_1 * B_1;
                            C_1_2 += A_1 * B_2;
                            C_2_0 += A_2 * B_0;
                            C_2_1 += A_2 * B_1;
                            C_2_2 += A_2 * B_2;
                        }
                        
                        C[i2 * k + j2] = C_0_0;
                        C[i2 * k + (j2 + 1)] = C_0_1;
                        C[i2 * k + (j2 + 2)] = C_0_2;
                        C[(i2 + 1) * k + j2] = C_1_0;
                        C[(i2 + 1) * k + (j2 + 1)] = C_1_1;
                        C[(i2 + 1) * k + (j2 + 2)] = C_1_2;
                        C[(i2 + 2) * k + j2] = C_2_0;
                        C[(i2 + 2) * k + (j2 + 1)] = C_2_1;
                        C[(i2 + 2) * k + (j2 + 2)] = C_2_2;
                    }
                }
            }
        }
    }
*/
    
    int i1 = i, j1 = j, k1 = k;
    int ni = i + b > n ? n : i + b;
    int nj = j + b > n ? n : j + b;
    int nk = k + b > n ? n : k + b;

    for (i1 = i; i1 < ni; i1 += 3)
    {
        for (j1 = j; j1 < nj; j1 += 3)
        {
            int t = i1 * n + j1;
            int tt = t + n;
            int ttt = tt + n;
            register double c00 = C[t];
            register double c01 = C[t + 1];
            register double c02 = C[t + 2];
            register double c10 = C[tt];
            register double c11 = C[tt + 1];
            register double c12 = C[tt + 2];
            register double c20 = C[ttt];
            register double c21 = C[ttt + 1];
            register double c22 = C[ttt + 2];

            for (k1 = k; k1 < nk; k1 += 3)
            {
		int l;
                for (l = 0; l < 3; l++)
                {
                    int ta = i1 * n + k1 + l;
                    int tta = ta + n;
                    int ttta = tta + n;
                    int tb = k1 * n + j1 + l * n;
                    register double a0 = A[ta];
                    register double a1 = A[tta];
                    register double a2 = A[ttta];
                    register double b0 = B[tb];
                    register double b1 = B[tb + 1];
                    register double b2 = B[tb + 2];

                    c00 -= a0 * b0;
                    c01 -= a0 * b1;
                    c02 -= a0 * b2;
                    c10 -= a1 * b0;
                    c11 -= a1 * b1;
                    c12 -= a1 * b2;
                    c20 -= a2 * b0;
                    c21 -= a2 * b1;
                    c22 -= a2 * b2;
                }
            }
            C[t] = c00;
            C[t + 1] = c01;
            C[t + 2] = c02;
            C[tt] = c10;
            C[tt + 1] = c11;
            C[tt + 2] = c12;
            C[ttt] = c20;
            C[ttt + 1] = c21;
            C[ttt + 2] = c22;
        }
    }
    return;
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    int ib, i, j, k, maxIndex;
    double max, sum;
    double *temprow = (double *) malloc(sizeof(double) * n);

    for (ib = 0; ib < n; ib += b)
    {
        for (i = ib; i < ib + b && i < n; i++)
        {
            // pivoting
            maxIndex = i;
            max = fabs(A[i * n + i]);
            
            for (j = i + 1; j < n; j++)
            {
                if (fabs(A[j * n + i]) > max)
                {
                    maxIndex = j;
                    max = fabs(A[j * n + i]);
                }
            }
            
            // if the matrix is singular
            if (max == 0)
            {
                return -1;
            }
            else
            {
                if (maxIndex != i)
                {
                    // save pivoting information
                    int temp = ipiv[i];
                    ipiv[i] = ipiv[maxIndex];
                    ipiv[maxIndex] = temp;
                    // swap rows
                    memcpy(temprow, A + i * n, n * sizeof(double));
                    memcpy(A + i * n, A + maxIndex * n, n * sizeof(double));
                    memcpy(A + maxIndex * n, temprow, n * sizeof(double));
                }
            }

            // factorization
            for (j = i + 1; j < n; j++)
            {
                A[j * n + i] = A[j * n + i] / A[i * n + i];
                
                for (k = i + 1; k < ib + b && k < n; k++)
                {
                    A[j * n + k] -= A[j * n + i] * A[i * n + k];
                }
            }
        }

        // update A(ib:end, end+1:n)
        for (i = ib; i < ib + b && i < n; i++)
        {
            for (j = ib + b; j < n; j++)
            {
                sum = 0;
                for (k = ib; k < i; k++)
                {
                    sum += A[i * n + k] * A[k * n + j];
                }
                A[i * n + j] -= sum;
            }
        }

        // update A(end+1:n, end+1:n)
        for (i = ib + b; i < n; i += b)
        {
            for (j = ib + b; j < n; j += b)
            {
                mydgemm(A, A, A, n, i, j, ib, b);
            }
        }
    }
    return 0;
}

