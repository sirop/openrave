/** \example orc_bodyfunctions.cpp
    \author Rosen Diankov

    Shows how to change body transforms, remove bodies, and get link information and change geometry properties.
 */
#include <openrave-core_c.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <sys/time.h>

#include <signal.h>

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

struct jointTime {
   double  joints[5];
   unsigned short  timeDiff;
};

void sigfunc(int sig) {
   int c;

   if(sig != SIGINT)
      return;
   else {
      printf("\n Stop the programm.\n");
   }
}


int main(int argc, char ** argv)
{
    ORCInitialize(0, 0); // no plugins, debuglevel=0
    void* env = ORCEnvironmentCreate();
    char* scenename= argv[1];
    ORCEnvironmentLoad(env, scenename);

    void* bodyhead = ORCEnvironmentGetKinBody(env, "HeadChassis");
    OpenRAVEReal values[5]; // = {0.0,PI/18, PI/4, 0.0, 0.0};
    OpenRAVEReal jacobian[15];
    double dchRollW = 0.0; // change to 0.5 ?
    double dchTiltW = 0.0;
    int m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info, rank;
    float rcond = -1.0;
    double ds[M];
    double db[LDB*N];
    double da[LDA*M];

    unsigned long long starttime, timediff;
    struct jointTime maxTime;
    maxTime.timeDiff = 0;

    signal(SIGINT, sigfunc);
   
    int chr, cht, p, r, t ; // index variables for the joints
    for ( chr = -40; chr <= 40; chr ++) {
        for ( cht = -40; cht <= 40; cht ++) {
            for ( p = -360; p <= 360; p ++)    {
                for ( r = -360; r <= 360; r ++)   {
                     for ( t = -360; t <= 360; t ++)  {
                          starttime = getMicrotime();

                          values[0] = (double)((chr)/2)*PI/180;
                          values[1] = (double)((cht-40)/2)*PI/180;
                          values[2] = (double)((p)/2)*PI/180;
                          values[3] = (double)((r)/2)*PI/180;
                          values[4] = (double)((t)/2)*PI/180;
                          ORCBodyJacobianAtDOFValues(bodyhead, values, jacobian );
                          
                          db[0] =  0.5 -jacobian[0]*dchRollW-jacobian[1]*dchTiltW;
                          db[1] =  0.5 -jacobian[5]*dchRollW-jacobian[6]*dchTiltW ;
                          db[2] =  0.5 -jacobian[10]*dchRollW-jacobian[11]*dchTiltW ;

                          memcpy(&da[0], &jacobian[2],  3*sizeof(double)); // , jacobian[3], jacobian[4],
                          memcpy(&da[3], &jacobian[7],  3*sizeof(double)); //, jacobian[8], jacobian[9],
                          memcpy(&da[6], &jacobian[12], 3*sizeof(double)); //jacobian[13], jacobian[14]
                           
                          info = LAPACKE_dgelsd( LAPACK_ROW_MAJOR, m, n, nrhs, da, lda, db, ldb,
                                 ds, rcond, &rank );
                          timediff = getMicrotime() - starttime;
                          //printf("LAPACKE_dgelsd time passed: %d\n",getMicrotime()-starttime);
                         /* Check for convergence */
                         if( info > 0 ) {
                              printf( "The algorithm computing SVD failed to converge;\n" );
                              printf( "the least squares solution could not be computed.\n" );
                         }
                         else {
                               printf("LAPACKE_dgelsd time passed: %d", (int)timediff);
                                /* Print minimum norm solution */
                               print_matrixd( "Minimum norm solution", n, nrhs, db, ldb );
                               if (timediff > maxTime.timeDiff) {
                                   maxTime.timeDiff = timediff;
                                   memcpy(maxTime.joints, values, 5*sizeof(double));
                               }
                         }
                     }
                }
            }
        }
    }
    
    
    int i;
    printf("Joint values (DEG) for maxTime %d: \n", maxTime.timeDiff);
    for (i=0; i<5; i++)
        printf("J[%d]= %f  ", i,  (maxTime.joints[i])/PI*180);
    printf("\n");
  
    ORCInterfaceRelease(bodyhead);
    ORCEnvironmentDestroy(env);
    ORCEnvironmentRelease(env);
    ORCDestroy();
    return 0;
};

/* Auxiliary routine: printing a matrix */
void print_matrixd( char const* desc, int m, int n, double* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
         printf("\n");
}

/**
 * Returns the current time in microseconds.
 */
long getMicrotime(){
	struct timeval currentTime;
	gettimeofday(&currentTime, NULL);
	return currentTime.tv_sec * (int)1e6 + currentTime.tv_usec;
}