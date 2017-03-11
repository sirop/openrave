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

#define M 4
#define N 3
#define NRHS 1
#define LDA N
#define LDB NRHS

//#define NumSteps 100

/* Auxiliary routines prototypes */
extern void print_matrixd( char const* desc, int m, int n, double* a, int lda );
extern long getMicrotime();
extern void  euler2quat(double ai, double aj,  double ak, double quaternion[4]);
extern void qmult(double q1[4], double q2[4], double qRes[4]);

static double fPanA = 0;
static double fTiltA = 0;
static double fRollA = 0;

static double fPanB = 0;
static double fTiltB = 0;
static double fRollB = 0;

static int deltaPanA = 0;
static int deltaTiltA = 0;
static int deltaRollA = 0;

static int deltaPanB = 0;
static int deltaTiltB = 0;
static int deltaRollB = 0;



static double panOffset = 0;
static double rollOffset = -0.2677953/180.0f*PI;
static double tiltOffset = 5.300267/180.0f*PI;



static double fperiod = 1.0/2000.0;

#define MaxRollVel (673474)  //  Max Roll command velocity in 0.1Motor counts /sec
#define MaxPanVel (673475)  //  Max Pan command velocity in 0.1Motor counts /sec
#define MaxTiltVel (133119)  //  Max Tilt command velocity in 0.1Motor counts /sec

#define RollVel (0)  //  Max Roll command velocity in 0.1Motor counts /sec
#define PanVel (30000)  //  Max Pan command velocity in 0.1Motor counts /sec
#define TiltVel (0)  //  Max Tilt command velocity in 0.1Motor counts /sec

#define StepScalePan  1.0
#define StepScaleRoll 1.0
#define StepScaleTilt 1.0

#define PAN_JOINTDEG2COUNTS  418.84444 //  #  counts/deg
#define ROLL_JOINTDEG2COUNTS 94.7911111
#define TILT_JOINTDEG2COUNTS 73.9555555556   // 26*1024/360 = 73.9555555556



int main(int argc, char ** argv)
{
    ORCInitialize(0, 1); // no plugins, debuglevel=1
    void* env = ORCEnvironmentCreate();
    char* scenename= argv[1];
    int NumSteps = atoi(argv[2]);
    ORCEnvironmentLoad(env, scenename);

    double pan = (double) (deltaPanA/PAN_JOINTDEG2COUNTS)/180*PI;
    double roll = (double)(deltaRollA/ROLL_JOINTDEG2COUNTS)/180*PI + rollOffset;
    double tilt = (double) (deltaTiltA/TILT_JOINTDEG2COUNTS)/180*PI + tiltOffset;


    void* bodyhead = ORCEnvironmentGetKinBody(env, "HeadChassis");
    OpenRAVEReal values[5] = { rollOffset, tiltOffset, pan, roll, tilt};
    OpenRAVEReal jacobian[20];
    
    int m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info, rank;
    float rcond = -1.0;
    int i = 0;
    // for (i=0; i<20; i++)
        // printf("J[i]= %f  ", jacobian[i]);
    // printf("\n");
    for ( i=0; i < NumSteps; i++) {
         unsigned long long starttime = getMicrotime();
         double dchRollW = 0.0;
         double dchTiltW = 0.0;
          
         //ORCBodyGetDOFValues(bodyhead, values);

 
        

         pan = values[2];
         roll = values[3];
         tilt = values[4];

         double quatEE [4] = {0, 0, 0, 0};
         double quatDelta [4] = {0, 0, 0, 0};
         double quatRes [4] = {0, 0, 0, 0};

         double StepPan = (double)  (PanVel)/MaxPanVel;  // velocity as a % of full speed
         double StepRoll = (double)  (RollVel)/MaxRollVel;  // velocity as a % of full speed
         double StepTilt = (double)  (TiltVel )/MaxTiltVel;  // velocity as a % of full speed
         StepPan = StepPan*PI*fperiod*StepScalePan;// Percentage
         StepRoll = StepRoll*PI*fperiod*StepScaleRoll; 
         StepTilt = StepTilt*PI*fperiod*StepScaleTilt;   
           
         euler2quat(StepRoll, StepTilt, StepPan, quatDelta);

         ORCBodyRotationJacobianAtDOFValues(bodyhead, values, jacobian );
         ORCBodyGetEEQuaternion(bodyhead, quatEE);
         qmult(quatEE, quatDelta, quatRes);

         double db[LDB*M] = { quatRes[0]-jacobian[0]*dchRollW-jacobian[1]*dchTiltW,
                              quatRes[1]-jacobian[5]*dchRollW-jacobian[6]*dchTiltW,
                              quatRes[2]-jacobian[10]*dchRollW-jacobian[11]*dchTiltW,
                              quatRes[3]-jacobian[15]*dchRollW-jacobian[16]*dchTiltW,
                        };
         double da[LDA*M] = { jacobian[2], jacobian[3], jacobian[4],
                              jacobian[7], jacobian[8], jacobian[9],
                              jacobian[12], jacobian[13], jacobian[14],
                              jacobian[17], jacobian[18], jacobian[19],
                        };
         double ds[M];
    
        info = LAPACKE_dgelsd( LAPACK_ROW_MAJOR, m, n, nrhs, da, lda, db, ldb,
                        ds, rcond, &rank );
        /*info = LAPACKE_dgels( LAPACK_ROW_MAJOR, 'N', m, n, nrhs, da, lda,
                        db, ldb );*/
         printf("LAPACKE_dgelsd time passed: %d\n",getMicrotime()-starttime);
       /* Check for convergence */
          if( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge;\n" );
                printf( "the least squares solution could not be computed.\n" );
                exit( 1 );
          }
        /* Print minimum norm solution */
          print_matrixd( "Minimum norm solution", n, nrhs, db, ldb );
          fPanA = (db[0]/PI*180.0f*PAN_JOINTDEG2COUNTS) + fPanA ; // fPanA is decimal accumulator
          fPanB = (db[0]/PI*180.0f*PAN_JOINTDEG2COUNTS) + fPanB ;
          //fStepPan = fPanA;
          deltaPanA = floor(fPanA) + deltaPanA;
          deltaPanB = floor(fPanB) + deltaPanB;
          fPanA = fPanA - floor(fPanA);
          fPanB = fPanB - floor(fPanB);            

          fRollA = (db[1]/PI*180.0f*ROLL_JOINTDEG2COUNTS) + fRollA;
          fRollB = (db[1]/PI*180.0f*ROLL_JOINTDEG2COUNTS) + fRollB;
          //fStepRoll = fRollA;
          deltaRollA = floor(fRollA) + deltaRollA ;
          deltaRollB = floor(fRollB) + deltaRollB ;
          fRollA = fRollA - floor(fRollA);
          fRollB = fRollB - floor(fRollB);
       
          fTiltA = (db[2]/PI*180.0f*TILT_JOINTDEG2COUNTS) + fTiltA + StepTilt  ;
          fTiltB = (db[2]/PI*180.0f*TILT_JOINTDEG2COUNTS) + fTiltB + StepTilt  ;
          //fStepTilt = fTiltA;
          deltaTiltA =  floor(fTiltA) + deltaTiltA  ;
          deltaTiltB =  floor(fTiltB) + deltaTiltB ;
          fTiltA = fTiltA - floor(fTiltA);
          fTiltB = fTiltB - floor(fTiltB);
          /*db[0] = floor(db[0]*1e+8)*1e-8;
          db[1] = floor(db[1]*1e+8)*1e-8;
          db[2] = floor(db[2]*1e+8)*1e-8;*/
          values[2] = db[0]+values[2]; //  to deltaTiltA ?
          values[3] = db[1]+values[3];
          values[4] = db[2]+values[4];
          printf("Pan: %10.22e, Roll %10.22e, Tilt: %10.22e  [DEG] \n", values[2]/PI*180.0f, values[3]/PI*180.0f, values[4]/PI*180.0f );
           printf("Pan: %d, Roll %d, Tilt: %d  [Counts] \n", deltaPanA, deltaRollA, deltaTiltA );
          ORCBodySetDOFValues( bodyhead, values);
    }

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
                for( j = 0; j < n; j++ ) printf( " %6.22e", a[i*lda+j] );
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

void  euler2quat(double ai, double aj,  double ak, double quaternion[4])
{
    /***
    Return `quaternion` from Euler angles and axis sequence `axes`
    ***/

    // axis sequences for Euler angles
    int _NEXT_AXIS[4] = {1, 2, 0, 1};
    int _AXES2TUPLE[4] = {0, 0, 0, 0}; // 'sxyz'
    int firstaxis = _AXES2TUPLE[0];
    int parity = _AXES2TUPLE[1];
    int repetition = _AXES2TUPLE[2];
    int frame = _AXES2TUPLE[3];

    int i = firstaxis + 1;
    int j = _NEXT_AXIS[i+parity-1] + 1;
    int k = _NEXT_AXIS[i-parity] + 1;

    if (frame){
        double temp = ak;
        ak = ai;
        ai = temp;
    }
    if (parity)
        aj = -aj;

    ai /= 2.0;
    aj /= 2.0;
    ak /= 2.0;
    double ci = cos(ai); // speedup with sincos
    double si = sin(ai);
    double cj = cos(aj);
    double sj = sin(aj);
    double ck = cos(ak);
    double sk = sin(ak);
    double cc = ci*ck;
    double cs = ci*sk;
    double sc = si*ck;
    double ss = si*sk;

    
    if (repetition) {
        quaternion[0] = cj*(cc - ss);
        quaternion[i] = cj*(cs + sc);
        quaternion[j] = sj*(cc + ss);
        quaternion[k] = sj*(cs - sc);
    }
    else{
        quaternion[0] = cj*cc + sj*ss;
        quaternion[i] = cj*sc - sj*cs;
        quaternion[j] = cj*ss + sj*cc;
        quaternion[k] = cj*cs - sj*sc;
    }
    if (parity)
        quaternion[j] *= -1.0;
}

void qmult(double q1[4], double q2[4], double qRes[4])
{
    double w1 = q1[0];
    double x1 = q1[1];
    double y1 = q1[2];
    double z1 = q1[3];

    double w2 = q2[0];
    double x2 = q2[1];
    double y2 = q2[2];
    double z2 = q2[3];

    qRes[0] = w1*w2 - x1*x2 - y1*y2 - z1*z2;
    qRes[1] = w1*x2 + x1*w2 + y1*z2 - z1*y2;
    qRes[2] = w1*y2 + y1*w2 + z1*x2 - x1*z2;
    qRes[3] = w1*z2 + z1*w2 + x1*y2 - y1*x2;
}