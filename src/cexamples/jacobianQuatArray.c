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
#include <time.h>

#include <signal.h>

//#include  "lapacke.h"

//#define PI  3.1415926535
#ifndef PI
#define PI 3.14159265358979323846
#endif

#define M 4
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
   //int c;

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
    double quatEE [4] = {0, 0, 0, 0};
    OpenRAVEReal jacobian[20];
    /* double dchRollW = 0.0; // change to 0.5 ?
    double dchTiltW = 0.0;
    int m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info, rank;
    float rcond = -1.0;
    double ds[M];
    double db[LDB*M];
    double da[LDA*M];
 */
   /*  unsigned long long starttime, timediff;
    struct jointTime maxTime;
    maxTime.timeDiff = 0; */

    signal(SIGINT, sigfunc);

    FILE *logfile;
    char logfilestring[22]="logfile";
    time_t now = time(NULL);
    struct tm *ti = localtime(&now);
    sprintf(&logfilestring[8], "%02d%02d%04d%02d%02d", ti->tm_mday, ti->tm_mon+1, 
    ti->tm_year + 1900, ti->tm_hour,ti->tm_min);
    logfile = fopen( logfilestring,  "w");

     unsigned long long counter = 0;
    
    int chr, cht, p, r, t ; // index variables for the joints
    for ( chr = -2; chr <= 2; chr ++) {
        for ( cht = -2; cht <= 2; cht ++) {
            for ( p = -91; p <= 91; p ++)    {
                for ( r = -4; r <= 4; r ++)   {
                     for ( t = -2; t <= 2; t ++)  {
                          //starttime = getMicrotime();

                          values[0] = (((double)chr)/2.0f)*PI/180.0f;
                          values[1] = (((double)cht)/2.0f)*PI/180.0f;
                          values[2] = (((double)p)/2.0f)*PI/180.0f;
                          values[3] = (((double)r)/2.0f)*PI/180.0f;
                          values[4] = (((double)t)/2.0f)*PI/180.0f;
                          
                          counter++;
                          if ((counter % 1000) == 0)
                              printf("%.22f %.22f %.22f %.22f %.22f counter = %llu\n",
                              (double)values[0], (double)values[1], (double)values[2],(double) values[3], (double)values[4], counter);

                          ORCBodyRotationJacobianAtDOFValues(bodyhead, values, jacobian );
                          ORCBodyGetEEQuaternion(bodyhead, quatEE);
                          fprintf(logfile, "%.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f \
                                            %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f %.22f \n",
                                            values[0], values[1], values[2], values[3], values[4],
                                            jacobian[0],  jacobian[1],  jacobian[2],  jacobian[3],  jacobian[4],
                                            jacobian[5],  jacobian[6],  jacobian[7],  jacobian[8],  jacobian[9],
                                            jacobian[10],  jacobian[11],  jacobian[12],  jacobian[13],  jacobian[14],
                                            jacobian[15],  jacobian[16],  jacobian[17],  jacobian[18],  jacobian[19],
                                            quatEE[0],  quatEE[1],  quatEE[2],  quatEE[3]
                           );
                          
                         /* db[0] =  0.5 -jacobian[0]*dchRollW-jacobian[1]*dchTiltW;
                          db[1] =  0.5 -jacobian[5]*dchRollW-jacobian[6]*dchTiltW ;
                          db[2] =  0.5 -jacobian[10]*dchRollW-jacobian[11]*dchTiltW ;

                          memcpy(&da[0], &jacobian[2],  3*sizeof(double)); // , jacobian[3], jacobian[4],
                          memcpy(&da[3], &jacobian[7],  3*sizeof(double)); //, jacobian[8], jacobian[9],
                          memcpy(&da[6], &jacobian[12], 3*sizeof(double)); //jacobian[13], jacobian[14]
                           
                          info = LAPACKE_dgelsd( LAPACK_ROW_MAJOR, m, n, nrhs, da, lda, db, ldb,
                                 ds, rcond, &rank );
                          timediff = getMicrotime() - starttime;*/
                          //printf("LAPACKE_dgelsd time passed: %d\n",getMicrotime()-starttime);
                         /* Check for convergence */
                        /* if( info > 0 ) {
                              printf( "The algorithm computing SVD failed to converge;\n" );
                              printf( "the least squares solution could not be computed.\n" );
                         }
                         else {
                               printf("LAPACKE_dgelsd time passed: %d", (int)timediff);
                               print_matrixd( "Minimum norm solution", n, nrhs, db, ldb );
                               if (timediff > maxTime.timeDiff) {
                                   maxTime.timeDiff = timediff;
                                   memcpy(maxTime.joints, values, 5*sizeof(double));
                               }*/
                         
                     }
                }
            }
        }
    }
    
    fclose(logfile);
    /* int i;
    printf("Joint values (DEG) for maxTime %d: \n", maxTime.timeDiff);
    for (i=0; i<5; i++)
        printf("J[%d]= %f  ", i,  (maxTime.joints[i])/PI*180);
    printf("\n"); */
  
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