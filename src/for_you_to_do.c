#include "../include/for_you_to_do.h"
#include <math.h>

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 128;
  
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
    double * tmp_row = (double *)malloc(sizeof(double) * n);
    
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
        
        // if the matrix A is singular
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
    double sum = 0;
    double * tmp_B = (double *)malloc(sizeof(double) * n);
    
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
    return 0;
}

