/*
	# File:    	qr.c
	# Purpose: 	Demonstrate the QR decomposition. The method of Householder
	|			reflections should be used:
	# Input:	rectangular MxN matrix A from the file A.txt which is the
	|			matrix of linear algebraic equasions system. The input file
	|			contains the M and N numbers at the first line. Below them there
	|			are the coefficients of the matrix.
	# Output:   the Q matrix, the R matrix and the QR multiplication (A matrix)
	# Compile:  Using gcc (can also be built it via makefile)
	|           gcc -g -Wall -o qr qr.c
	|           OR make all
	# Usage:    ./qr
	# Note:		The input matrix is hardcoded now. If you want to enter your
	|			own MxN matrix, write the corresponding function to fill A
	|			before calculations
	# Author: 	Nikolai Gaiduchenko, MIPT 2017, 513 group
*/
//==============================================================================
// INCLUDE SECTION
//==============================================================================

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <c6x.h>
#include <ti/csl/csl_chip.h>
#include <ti/csl/csl_tsc.h>
#include <ti/csl/csl_cacheAux.h>
//#include <DSPF_sp_mat_mul\c66\DSPF_sp_mat_mul.h>
//#include <DSPF_sp_mat_trans\c66\DSPF_sp_mat_trans.h>
#include "DSPF_sp_mat_mul.h"
#include "DSPF_sp_mat_trans.h"
#include "DSPF_sp_vecsum_sq.h"
#include "DSPF_sp_vecadd.h"
#include "Timer.h"
#include "InitPLL.h"

extern int * ResultAddr; // 存放计数的地址
float result;	// 存放计数结果
int Fre = 1000;	// 频率(MHz)
/*
unsigned char DSP_NUM;
DSP_NUM =( unsigned char)((*(unsigned int volatile *)(0x02320020)) & 0x000E);
DSP_NUM = DSP_NUM >> 1;*/

//==============================================================================
// MATRIX (STRUCT AND INPUT) SECTION
//==============================================================================

typedef struct matrix_type {
	int m, n;
	float ** v;
} mat_t, *mat;

#define M 64
#define N 32
float Mat[M][N];
mat qgs_Q;

void SysInit();
//==============================================================================
// MATRIX FUNCTIONS PROTOTYPES
//==============================================================================

mat matrix_new(int m, int n);
// PURPOSE:			allocate memory for the matrix
// RETURN VALUE:	the matrix created

void matrix_delete(mat m);
// PURPOSE:			free the memory being allocated to store the matrix
// RETURN VALUE:	Nothing

void matrix_transpose(mat m);
// PURPOSE:			Transpose the matrix (A[i][j] = A[j][i])
// RETURN VALUE:	Nothing

mat matrix_copy(int n, float **a, int m);
// PURPOSE:			Copy the matrix given w/ memory allocation
// RETURN VALUE:	The new matrix which is a copy of the matrix given

mat matrix_mul(mat x, mat y);
// PURPOSE:			Multiply the matrices given
// RETURN VALUE:	The new matrix as a result of multiplication

mat matrix_minor(mat x, int d);
// PURPOSE:			Calculate the Dth minor of the matrix
// RETURN VALUE:	The new matrix whis is the Dth minor of the original one

void matrix_show(mat m);
// PURPOSE:			Print the matrix
// RETURN VALUE:	Nothing

//==============================================================================
// VECTOR FUNCTION PROTOTYPES
//==============================================================================

float *vmadd(float a[], float b[], float s, float c[], int n);
// RETURN VALUE:	c = a + b*s

float* mcol(mat m, float *v, int c);
// PURPOSE:			Take the c-th column from the matrix m and put it into
//					vector v
// RETURN VALUE:	Vector v which is the column c in matrix m

mat vmul(float v[], int n);
// PURPOSE:			Create matrix m = I - v * v^T,
//					where I is the identical matrix
// RETURN VALUE:	The matrix created

float vnorm(float x[], int n);
// PURPOSE:			Find the vector norm (module)
// RETURN VALUE:	||x||

float* vdiv(float x[], float d, float y[], int n);
// PURPOSE:			Divide the vector x[n] by d and store the result in y[n]
// RETURN VALUE:	y[]
void* aligned_malloc(size_t required_bytes, size_t alignment)

{

    int offset = alignment - 1 + sizeof(void*);

    void* p1 = (void*)malloc(required_bytes + offset);

    if (p1 == NULL)

        return NULL;

    void** p2 = (void**)( ( (size_t)p1 + offset ) & ~(alignment - 1) );

    p2[-1] = p1;

    return p2;

}

void aligned_free(void *p2)

{

    void* p1 = ((void**)p2)[-1];

    free(p1);

}
//==============================================================================
// THE HOUSEHOLDER METHOD IMPLEMENTATION
//==============================================================================

void householder(mat m, mat *R, mat *Q)
{
	mat *q = (mat*)aligned_malloc(m->m * sizeof(mat),128);
	mat z = m, z1;

	int k, i;
	#pragma MUST_ITERATE(1,32,1)
	for (k = 0; k < m->n && k < m->m - 1; k++) {

		float *e = (float*)malloc(m->m * sizeof(float));
		float *x = (float*)aligned_malloc(m->m * sizeof(float), 128);
		float a;
		z1 = matrix_minor(z, k);
		if (z != m) matrix_delete(z);
		z = z1;

		mcol(z, x, k);   //列向量x
		a = vnorm(x, m->m);
		if (m->v[k][k] > 0) a = -a;

		memset(e, 0, m->m * sizeof(float));
		e[k] = 1;
//		for (i = 0; i < m->m; i++)
//			e[i] = (i == k) ? 1 : 0;

		vmadd(x, e, a, e, m->m);
		vdiv(e, vnorm(e, m->m), e, m->m);
		q[k] = vmul(e, m->m);

		z1 = matrix_mul(q[k], z);

		if (z != m) matrix_delete(z);
		z = z1;

		free(e);
		aligned_free(x);
	}

	matrix_delete(z);
	*Q = q[0];
	*R = matrix_mul(q[0], m);

	#pragma MUST_ITERATE(1,,)
	for (i = 1; i < m->n && i < m->m - 1; i++) {
		z1 = matrix_mul(q[i], *Q);
		if (i > 1) matrix_delete(*Q);
		*Q = z1;
		matrix_delete(q[i]);
	}

	matrix_delete(q[0]);

	z = matrix_mul(*Q, m);

	matrix_delete(*R);
	*R = z;

	qgs_Q = matrix_new((*Q)->m, (*Q)->n);
	matrix_transpose(*Q);
}

//==============================================================================
// MATRIX FUNCTIONS IMPLEMENTATIONS
//==============================================================================
//dobule g_matrix_data = {1,2,3,4}
//mat_t g_matInput = {4,5,
mat matrix_new(int m, int n)
{
	int i;
	mat x = (mat)aligned_malloc(sizeof(mat_t), 128);
	x->v = (float **)aligned_malloc(sizeof(float) * m, 128);
	x->v[0] = (float *)aligned_malloc(sizeof(float) * m * n, 128);

	#pragma MUST_ITERATE(1,,);
	for (i = 0; i < m; i++)
		x->v[i] = x->v[0] + n * i;
	x->m = m;
	x->n = n;
	return x;
}

void matrix_delete(mat m)
{
	aligned_free(m->v[0]);
	aligned_free(m->v);
	aligned_free(m);
}

void matrix_transpose(mat m)
{
#if 1
	DSPF_sp_mat_trans(m->v[0], m->m, m->n, qgs_Q->v[0]);
#else
	int i, j;
	for (i = 0; i < m->m; i++) {
		for (j = 0; j < i; j++) {
			float t = m->v[i][j];
			m->v[i][j] = m->v[j][i];
			m->v[j][i] = t;
		}
	}
#endif
}

mat matrix_copy(int n, float **a, int m)
{
	int i, j;
	mat x = matrix_new(m, n);
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			x->v[i][j] = a[i][j];
	return x;
}


typedef struct {
    int job;
    int am;
    int an;
    int bn;
    float *a;
    float *b;
    float *c;
    int res;
}qr_task;

#pragma DATA_SECTION (g_task, ".global_shared");
qr_task g_task[8];

mat matrix_mul(mat x, mat y)
{
#if 1
/*
	mat r;

	r = matrix_new(x->m, y->n);

	DSPF_sp_mat_mul(x->v[0], x->m, x->n, y->v[0] ,y->n, (float *restrict)r->v[0]);*/
	int i;
    mat r;
    r = matrix_new(x->m, y->n);

    memset(g_task, 0, sizeof(qr_task)*8);

    g_task[0].b = g_task[1].b = g_task[2].b = g_task[3].b = y->v[0];
    g_task[0].bn = g_task[1].bn = g_task[2].bn = g_task[3].bn = y->n;
    g_task[0].an = g_task[1].an = g_task[2].an = g_task[3].an = x->n;
    g_task[0].am = g_task[1].am = g_task[2].am = g_task[3].am = x->m/4;

	#pragma MUST_ITERATE(2,,);
	#pragma UNROLL(2);
    for(i = 0; i < 4; i++)
    {
    	g_task[i].a = x->v[0] + i * x->m/4 * x->n;
    	g_task[i].c = r->v[0] + i * x->m/4 * y->n;
    }
/*    g_task[0].a = x->v[0];
    g_task[0].c = r->v[0];

	g_task[1].a = x->v[0] + g_task[0].am * g_task[0].an;
	g_task[1].c = r->v[0] + g_task[0].am * g_task[0].bn;

	g_task[2].a = x->v[0] + g_task[0].am * g_task[0].an + g_task[1].am * g_task[1].an;
	g_task[2].c = r->v[0] + g_task[0].am * g_task[0].bn + g_task[1].am * g_task[1].bn;

    g_task[3].a = x->v[0] + g_task[0].am * g_task[0].an + g_task[1].am * g_task[1].an + g_task[2].am * g_task[2].an;
    g_task[3].c = r->v[0] + g_task[0].am * g_task[0].bn + g_task[1].am * g_task[1].bn + g_task[2].am * g_task[2].bn;*/

    g_task[0].job = g_task[1].job = g_task[2].job = g_task[3].job = 1;

    CACHE_wbL1d ((void *)g_task[1].a, sizeof(float) * 8 * g_task[1].am * g_task[1].an, CACHE_WAIT);
    CACHE_wbL1d ((void *)g_task[1].b, sizeof(float) * 8 * g_task[1].an * g_task[1].bn, CACHE_WAIT);

    CACHE_wbL1d ((void *)g_task[2].a, sizeof(float) * 8 * g_task[2].am * g_task[2].an, CACHE_WAIT);
    CACHE_wbL1d ((void *)g_task[2].b, sizeof(float) * 8 * g_task[2].an * g_task[2].bn, CACHE_WAIT);

    CACHE_wbL1d ((void *)g_task[3].a, sizeof(float) * 8 * g_task[3].am * g_task[3].an, CACHE_WAIT);
    CACHE_wbL1d ((void *)g_task[3].b, sizeof(float) * 8 * g_task[3].an * g_task[3].bn, CACHE_WAIT);

    CACHE_wbL1d ((void *)&(g_task[1].job), sizeof(int), CACHE_WAIT);
    CACHE_wbL1d ((void *)&(g_task[2].job), sizeof(int), CACHE_WAIT);
    CACHE_wbL1d ((void *)&(g_task[3].job), sizeof(int), CACHE_WAIT);

    DSPF_sp_mat_mul(g_task[0].a, g_task[0].am, g_task[0].an, g_task[0].b, g_task[0].bn, g_task[0].c);

    while(_lor(_lor(g_task[1].job,g_task[2].job), g_task[3].job))
    {
       CACHE_invL1d ((void *)&(g_task[1].job), sizeof(int), CACHE_WAIT);
       CACHE_invL1d ((void *)&(g_task[2].job), sizeof(int), CACHE_WAIT);
       CACHE_invL1d ((void *)&(g_task[3].job), sizeof(int), CACHE_WAIT);
    }
    CACHE_invL1d ((void *)g_task[1].c, sizeof(float) * g_task[1].am * g_task[1].bn, CACHE_WAIT);
    CACHE_invL1d ((void *)g_task[2].c, sizeof(float) * g_task[2].am * g_task[2].bn, CACHE_WAIT);
    CACHE_invL1d ((void *)g_task[3].c, sizeof(float) * g_task[3].am * g_task[3].bn, CACHE_WAIT);

#else
	int i, j, k;
	mat r;
	if (x->n != y->m) return 0;
	r = matrix_new(x->m, y->n);
	for (i = 0; i < x->m; i++)
		for (j = 0; j < y->n; j++)
			for (k = 0; k < x->n; k++)
				r->v[i][j] += x->v[i][k] * y->v[k][j];
//				r->v[i][j] = 1;
#endif
	return r;
}

mat matrix_minor(mat x, int d)
{
	int i, j;
	mat m = matrix_new(x->m, x->n);
	for (i = 0; i < d; i++)
		m->v[i][i] = 1;
	#pragma MUST_ITERATE(2,,);
	#pragma UNROLL(2);
	for (i = d; i < x->m; i++)
		#pragma MUST_ITERATE(1,,)
		for (j = d; j < x->n; j += 4)
		{
			_amemd8(&(m->v[i][j])) = _amemd8_const(&(x->v[i][j]));
			_amemd8(&(m->v[i][j+2])) = _amemd8_const(&(x->v[i][j+2]));
		}
	return m;

}

//==============================================================================
// VECTOR FUNCTIONS IMPLEMENTATIONS
//==============================================================================

/* c = a + b * s */
float *vmadd(float a[], float b[], float s, float c[], int n)
{
#if 1
/*	int i;
	double ss = _ftod(s, s);
    _nassert(((int)a & 7) ==0);
    _nassert(((int)b & 7) ==0);
    _nassert(((int)c & 7) ==0);
    #pragma MUST_ITERATE(2,,2);
	#pragma UNROLL(2);
    for(i=0; i<n; i+=2) {
    	_amemd8(&c[i]) = _daddsp(_amemd8_const(&a[i]), _dmpysp(ss, _amemd8_const(&b[i])));
    }*/
    int    i;
    double x11x10, x21x20, x13x12, x23x22, m_m;

    _nassert(n > 0);
    _nassert(n % 4 == 0);
    _nassert((int)a % 8 == 0);
    _nassert((int)b % 8 == 0);
    _nassert((int)c  % 8 == 0);

    m_m = _ftod(s, s);
    #pragma MUST_ITERATE(1,,1)
    for (i = 0; i < n; i += 4)
    {
        x11x10 = _amemd8((void *) (a + i));
        x13x12 = _amemd8((void *) (a + i + 2));
        x21x20 = _amemd8((void *) (b + i));
        x23x22 = _amemd8((void *) (b + i + 2));

        _amemd8((void*)&c[i])   = _daddsp(x21x20, _dmpysp(m_m, x11x10));
        _amemd8((void*)&c[i+2]) = _daddsp(x23x22, _dmpysp(m_m, x13x12));
    }
	return c;
#else
	int i;
	_nassert(((int)a & 7) ==0);
	_nassert(((int)b & 7) ==0);
	_nassert(((int)c & 7) ==0);
	for (i = 0; i < n; i++)
		c[i] = a[i] + s * b[i];
	return c;
#endif
}

/* m = I - v v^T */
mat vmul(float v[], int n)
{
	int i, j;
	double v_v;
	float t = -2;
	double tt;
	tt = _ftod(t, t);
	mat x = matrix_new(n, n);
#if 0
	#pragma MUST_ITERATE(1,,)
	for (i = 0; i < n; i++)
	{
		#pragma MUST_ITERATE(1,,)
		for (j = 0; j < i; j++)
		{
			x->v[i][j] = x->v[j][i] = -2 * v[i] * v[j];
		}
		x->v[i][i] = -2 * v[i] * v[i] + 1;
	}
#else
	for (i = 0; i < n; i++)
	{
		v_v = _ftod(v[i], v[i]);
		for (j = 0; j < n; j += 2)
		{
			_amemd8(&(x->v[i][j])) = _dmpysp(tt,_dmpysp(_amemd8(&v[j]),v_v));
		}
		x->v[i][i]++;
	}
#endif
	return x;
}

/* ||x|| */
float vnorm(float x[], int n)
{
#if 0
	int i;
	float sum = 0;
	for (i = 0; i < n; i++) sum += x[i] * x[i];
	return sqrt(sum);
#else
	return sqrt(DSPF_sp_vecsum_sq(x, n));
#endif
}

/* y = x / d */
float* vdiv(float x[], float d, float y[], int n)
{
	int i;
	double dd = _ftod(1/d, 1/d);
	#pragma MUST_ITERATE(1,,);
	for (i = 0; i < n; i += 4)
	{
		_amemd8(&y[i]) = _dmpysp(_amemd8_const(&x[i]), dd);
		_amemd8(&y[i+2]) = _dmpysp(_amemd8_const(&x[i+2]), dd);
	}
	return y;
}

/* take c-th column of m, put in v */
float* mcol(mat m, float *v, int c)
{
	int i;
#if 0
	qgs_Q = matrix_new(m->m, m->n);
	matrix_transpose(m);
	for (i = 0; i < m->m; i++)
		v[i] = qgs_Q->v[c][i];
	matrix_delete(qgs_Q);
#else
	_nassert((int)v % 8 == 0 );

	#pragma MUST_ITERATE(1,,);
	for (i = 0; i < m->m; i += 4)
	{
		_amemd8(&v[i]) = _amemd8_const(&Mat[c][i]);
		_amemd8(&v[i+2]) = _amemd8_const(&Mat[c][i+2]);
	}
//		v[i] = m->v[i][c];
#endif
	return v;
}

void matrix_show(mat m)
{
	int i, j;
	for(i = 0; i < m->m; i++) {
		for (j = 0; j < m->n; j++) {
			printf(" %8.3f", m->v[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

//==============================================================================
// MAIN CODE SECTION
//==============================================================================
unsigned int DSP_NUM;

int main()
{
	// read the input MxN matrix from the file
	// M - num of strings
	// N - num of columns

/*	*(unsigned int volatile *)(0x02320010) |= 0x000e;      //定义为输入
	DSP_NUM =( unsigned char)((*(unsigned int volatile *)(0x02320020)) & 0x000E);
	DSP_NUM = DSP_NUM >> 1;*/

	DSP_NUM = CSL_chipReadReg (CSL_CHIP_DNUM);
	if(DSP_NUM == 0)
	{
		unsigned i, k;
		mat inmat;
		mat R, Q;
		mat m, C;
/*
		uint32_t marIndex;
	    for(marIndex = 0x80000000/16/1024/1024; marIndex < 0x90000000/16/1024/1024; marIndex++)
	    {
	        CACHE_disableCaching ((Uint8)marIndex);
	    }*/

		//输入前做一次转置
		for(i=0;i<N;i++)
		{
			for(k=0;k<M;k++)
			{
				Mat[i][k] = (rand()/(float)(RAND_MAX/10));
			}
		}

		SysInit();
		TIMER_ClkOpen();
		// Allocate memory for input matrix (dynamic) MxN
		inmat = matrix_new(M, N);
		//	TIMER_END();
		//	result = *( (int*) ResultAddr );
		// read matrix coefficients from the input file
		for (i = 0; i < M; i++)
		{
			for (k = 0; k < N; k++)
			{
				inmat->v[i][k] = Mat[k][i];     //转置
			}
		}
		puts("inmat"); matrix_show(inmat);
		TIMER_START();
		householder(inmat, &R, &Q);
		TIMER_END();
		result = *( (int*) ResultAddr );

		puts("Q"); matrix_show(qgs_Q);
		puts("R"); matrix_show(R);

		// to show their product is the input matrix
		m = matrix_mul(qgs_Q, R);
		puts("Q * R"); matrix_show(m);

		C = m;
		vmadd(m->v[0], inmat->v[0], -1, C->v[0], M*N);
		puts("Verification"); matrix_show(C);

		printf("time=%f μs\n", (float)result*6/Fre);

		matrix_delete(R);
		matrix_delete(Q);
		matrix_delete(qgs_Q);
		matrix_delete(m);
		matrix_delete(inmat);
		return 0;
	}
	else if(DSP_NUM == 1)
	{
		while(1)
		{
            CACHE_invL1d ((void *)&(g_task[1].job), sizeof(int), CACHE_WAIT);
            CACHE_invL1d ((void *)g_task[1].a, sizeof(float)*8*g_task[1].am*g_task[1].an, CACHE_WAIT);
            CACHE_invL1d ((void *)g_task[1].b, sizeof(float)*8*g_task[1].an*g_task[1].bn, CACHE_WAIT);
			if(g_task[1].job == 1)
			{
				DSPF_sp_mat_mul(g_task[1].a, g_task[1].am, g_task[1].an, g_task[1].b, g_task[1].bn, g_task[1].c);
				g_task[1].job = 0;
			}
            CACHE_wbL1d ((void *)&(g_task[1].job), sizeof(int), CACHE_WAIT);
            CACHE_wbL1d ((void *)g_task[1].c, sizeof(float) * g_task[1].am * g_task[1].bn, CACHE_WAIT);
		}
	}
	else if(DSP_NUM == 2)
	{
		while(1)
		{
            CACHE_invL1d ((void *)&(g_task[2].job), sizeof(int), CACHE_WAIT);
            CACHE_invL1d ((void *)g_task[2].a, sizeof(float)*8*g_task[2].am*g_task[2].an, CACHE_WAIT);
            CACHE_invL1d ((void *)g_task[2].b, sizeof(float)*8*g_task[2].an*g_task[2].bn, CACHE_WAIT);
			if(g_task[2].job == 1)
			{
				DSPF_sp_mat_mul(g_task[2].a, g_task[2].am, g_task[2].an, g_task[2].b, g_task[2].bn, g_task[2].c);
				g_task[2].job = 0;
			}
            CACHE_wbL1d ((void *)&(g_task[2].job), sizeof(int), CACHE_WAIT);
            CACHE_wbL1d ((void *)g_task[2].c, sizeof(float) * g_task[2].am * g_task[2].bn, CACHE_WAIT);
		}
	}
	else if(DSP_NUM == 3)
	{
		while(1)
		{
            CACHE_invL1d ((void *)&(g_task[3].job), sizeof(int), CACHE_WAIT);
            CACHE_invL1d ((void *)g_task[3].a, sizeof(float)*8*g_task[3].am*g_task[3].an, CACHE_WAIT);
            CACHE_invL1d ((void *)g_task[3].b, sizeof(float)*8*g_task[3].an*g_task[3].bn, CACHE_WAIT);
			if(g_task[3].job == 1)
			{
				DSPF_sp_mat_mul(g_task[3].a, g_task[3].am, g_task[3].an, g_task[3].b, g_task[3].bn, g_task[3].c);
				g_task[3].job = 0;
			}
            CACHE_wbL1d ((void *)&(g_task[3].job), sizeof(int), CACHE_WAIT);
            CACHE_wbL1d ((void *)g_task[3].c, sizeof(float) * g_task[3].am * g_task[3].bn, CACHE_WAIT);
		}
	}
}

void SysInit(){
	unsigned int i;
	//disable cache
	*(unsigned int volatile *)(0x01845044) = 1;
	do{
	i = (*(unsigned int volatile *)(0x01845044)&0x1);
	}while(i!=0);
	*(unsigned int volatile *)(0x01840040) = 0;

	printf("Test begin\n");
	printf("main PLL 1GHz\n");
	MainPLL(40,1,1,1);//1GHz

}
//==============================================================================
// USAGE EXAMPLE
//==============================================================================
/*
	Input file A.txt contains:
	5 3
	12.000  -51.000    4.000
	 6.000  167.000  -68.000
	-4.000   24.000  -41.000
	-1.000    1.000   -0.000
	 2.000   -0.000    3.000
    Output:
    Q:
    0.846   -0.391    0.343    0.082    0.078
    0.423    0.904   -0.029    0.026    0.045
   -0.282    0.170    0.933   -0.047   -0.137
   -0.071    0.014   -0.001    0.980   -0.184
    0.141   -0.017   -0.106   -0.171   -0.969
    R:
   14.177   20.667  -13.402
   -0.000  175.043  -70.080
    0.000    0.000  -35.202
   -0.000   -0.000   -0.000
    0.000    0.000   -0.000
    Q * R:
   12.000  -51.000    4.000
    6.000  167.000  -68.000
   -4.000   24.000  -41.000
   -1.000    1.000   -0.000
    2.000   -0.000    3.000
*/
