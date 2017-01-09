/** \example orikfilter.cpp
    \author Rosen Diankov

    Shows how to use set a custom inverse kinematics filter to add extra constraints.

    <b>Full Example Code:</b>
 */
#include <openrave-core.h>
#include <openrave/utils.h>
#include <vector>
#include <sstream>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

#include  "lapacke.h"

using namespace OpenRAVE;
using namespace std;

/* Auxiliary routines prototypes */
extern void print_matrix( char const* desc, int m, int n, float* a, int lda );


int main(int argc, char ** argv)
{
    if( argc < 2 ) {
        RAVELOG_INFO("Robot XML Path needed.\n");
        return 1;
    }
    string robotname = argv[1];
    //string scenefilename = "data/pr2test1.env.xml";
    RaveInitialize(true);
    EnvironmentBasePtr penv = RaveCreateEnvironment();
    penv->Load(robotname);

    vector<RobotBasePtr> vrobots;
    penv->GetRobots(vrobots);
    RobotBasePtr probot = vrobots.at(0);
    //probot->SetActiveManipulator("leftarm_torso");
    vector<RobotBase::ManipulatorPtr> vmanips;
    vmanips = probot->GetManipulators();
    //stringstream ssout;
    cout << "vmanips.size()=" << vmanips.size() <<   '\n' ;
    RobotBase::ManipulatorPtr pmanip = vmanips.at(0);

    //probot->SetActiveDOFs(pmanip->GetArmIndices());
    vector<dReal> jacobian;
    pmanip->CalculateAngularVelocityJacobian(jacobian);
   // stringstream ssout;
    cout << "jacobian.size()=" << jacobian.size() <<   '\n';
    for(int i=0; i<jacobian.size(); ++i)
        cout << jacobian.at(i) << '\n';    
    


#define M 3
#define N 3
#define NRHS 1
#define LDA N
#define LDB NRHS

    int m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info, rank;
    float rcond = -1.0;
    float b[LDB*N] = { 0.0f, 0.0f, 0.5f};
    float a[LDA*M] = { jacobian.at(2), jacobian.at(3), jacobian.at(4),
                       jacobian.at(7), jacobian.at(8), jacobian.at(9),
                       jacobian.at(12), jacobian.at(13), jacobian.at(14)
                     };
    float s[M];
    
    info = LAPACKE_sgelsd( LAPACK_ROW_MAJOR, m, n, nrhs, a, lda, b, ldb,
                        s, rcond, &rank );
     
     /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge;\n" );
                printf( "the least squares solution could not be computed.\n" );
                exit( 1 );
        }
        /* Print minimum norm solution */
        print_matrix( "Minimum norm solution", n, nrhs, b, ldb );
        /* Print effective rank */
        printf( "\n Effective rank = %6i\n", rank );
        /* Print singular values */
        print_matrix( "Singular values", 1, m, s, 1 );
        //exit( 0 );    



    RaveDestroy();
    return 0;
}

/* Auxiliary routine: printing a matrix */
void print_matrix( char const* desc, int m, int n, float* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}
