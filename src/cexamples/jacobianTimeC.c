/** \example orc_bodyfunctions.cpp
    \author Rosen Diankov

    Shows how to change body transforms, remove bodies, and get link information and change geometry properties.
 */
#include <openrave-core_c.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <sys/time.h>

#include  "lapacke.h"

#define PI  3.1415926535

#define M 3
#define N 3
#define NRHS 1
#define LDA N
#define LDB NRHS

/* Auxiliary routines prototypes */
extern void print_matrixd( char const* desc, int m, int n, double* a, int lda );
extern long getMicrotime();

int main(int argc, char ** argv)
{
    ORCInitialize(0, 1); // no plugins, debuglevel=1
    void* env = ORCEnvironmentCreate();
    char* scenename= argv[1];
    ORCEnvironmentLoad(env, scenename);
    //ORCEnvironmentSetViewer(env,"qtcoin");

   // ORCSetDebugLevel(4); // set to debug

    //bool pole2added = true;
    void* bodyhead = ORCEnvironmentGetKinBody(env, "HeadChassis");
    OpenRAVEReal values[5] = {0.0,PI/18, PI/4, 0.0, 0.0};
    OpenRAVEReal jacobian[15];
    ORCBodyJacobianAtDOFValues(bodyhead, values, jacobian );
    int i;
    for (i=0; i<15; i++)
        printf("J[i]= %f  ", jacobian[i]);
    printf("\n");
    unsigned long long starttime = getMicrotime();
     double dchRollW = 0.5;
     double dchTiltW = 0.5;
     int m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info, rank;
     float rcond = -1.0;

     double db[LDB*N] = { 0.0-jacobian[0]*dchRollW-jacobian[1]*dchTiltW,
                             0.0-jacobian[5]*dchRollW-jacobian[6]*dchTiltW,
                             0.5-jacobian[10]*dchRollW-jacobian[11]*dchTiltW
                        };
     double da[LDA*M] = { jacobian[2], jacobian[3], jacobian[4],
                          jacobian[7], jacobian[8], jacobian[9],
                          jacobian[12], jacobian[13], jacobian[14]
                        };
       double ds[M];
    
       info = LAPACKE_dgelsd( LAPACK_ROW_MAJOR, m, n, nrhs, da, lda, db, ldb,
                        ds, rcond, &rank );
       printf("LAPACKE_dgelsd time passed: %d\n",getMicrotime()-starttime);
       /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge;\n" );
                printf( "the least squares solution could not be computed.\n" );
                exit( 1 );
        }
        /* Print minimum norm solution */
        print_matrixd( "Minimum norm solution", n, nrhs, db, ldb );
  
    ORCInterfaceRelease(bodyhead);
    ORCEnvironmentDestroy(env);
    ORCEnvironmentRelease(env);
    ORCDestroy();
};

/* Auxiliary routine: printing a matrix */
void print_matrixd( char const* desc, int m, int n, double* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}

/**
 * Returns the current time in microseconds.
 */
long getMicrotime(){
	struct timeval currentTime;
	gettimeofday(&currentTime, NULL);
	return currentTime.tv_sec * (int)1e6 + currentTime.tv_usec;
}