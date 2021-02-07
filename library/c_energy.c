#include "mex.h"
#include<math.h>

void kenergy(double *ek, double* At, const double * restrict vx, const double * restrict vy, const double * restrict vz, 
             const int ns, const double * restrict np, const double * restrict mass, const double cv, double bx0, double by0, double bz0){

  
  double ek_tmp, vs, vspara;
  const double csq = cv*cv;
  const double b0 = sqrt(bx0*bx0+by0*by0+bz0*bz0);

  int ifirst=0;
  int ilast=0;  

  for(int is=0;is<ns;is++){
    ifirst = ilast;
    ilast = ifirst + (int64_T) np[is];

    ek_tmp = 0.0;
    vs = 0.0;
    vspara = 0.0;

    for(int ip=ifirst;ip<ilast;ip++){
      ek_tmp += (cv/sqrt(csq-vx[ip]*vx[ip]-vy[ip]*vy[ip]-vz[ip]*vz[ip])-1.0)*mass[is];
      vs += vx[ip]*vx[ip] + vy[ip]*vy[ip]+vz[ip]*vz[ip];
      vspara += pow((vx[ip]*bx0 + vy[ip]*by0 + vz[ip]*bz0)/b0, 2);
    }  
     ek[is] = ek_tmp;
      //Temperature Anisotropy  T_perp/T_para
     At[is] = 0.5*(vs/vspara -1.0);        
  }
}

void mexFunction( int nlhs, mxArray *plhs[],       // �߂�l�̌�, �߂�l
                  int nrhs, const mxArray *prhs[]) // �����̌�, ����
{
     /* check for proper number of arguments */
    if(nrhs!=10) {
        mexErrMsgIdAndTxt("c_current:RHS","10 inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("c_current:LHS","3 outputs required.");
    }
    
    /* make sure the first 3 argument is 1d double vector */
    for(int i=0;i<3;i++){
      if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxGetNumberOfDimensions(prhs[i])!=2 ) {
          mexErrMsgIdAndTxt("c_current:particle", "Input multiplier must be a 1d double vector.");
      }
    }

    /* make sure the 4th input argument is scalar */
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3])!=1 ) {
        mexErrMsgIdAndTxt("c_current:ns","Input multiplier must be a scalar.");
    }

    /* make sure the 5th input argument 1d int vector */
    if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetNumberOfDimensions(prhs[4])!=2 ) {
        mexErrMsgIdAndTxt("c_current:np", "Input multiplier must be a 1d integer vector.");
    }

    /* make sure the 6th input argument is vector */
    if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mxGetNumberOfDimensions(prhs[5])!=2 ) {
        mexErrMsgIdAndTxt("c_current:mass","Input multiplier must be a 1d vector.");
    }

    /* make sure the 7th-10th input argument is scalar */
    for(int i=6;i<10;i++){
      if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxGetNumberOfElements(prhs[i])!=1 ) {
          mexErrMsgIdAndTxt("c_current:field size","Input multiplier must be a scalar.");
      }
    }
    /*�����̍s��ƃT�C�Y���擾����*/    
#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *vx, *vy, *vz;
    int ns;
    mxDouble *np, *mass;
    double cv, bx0, by0, bz0;
    
    vx = mxGetDoubles(prhs[0]);
    vy = mxGetDoubles(prhs[1]);
    vz = mxGetDoubles(prhs[2]);
    ns = mxGetScalar(prhs[3]);
    np = mxGetDoubles(prhs[4]);
    mass = mxGetDoubles(prhs[5]);
    cv = mxGetScalar(prhs[6]);
    bx0 = mxGetScalar(prhs[7]);
    by0 = mxGetScalar(prhs[8]);
    bz0 = mxGetScalar(prhs[9]);
#else
    double *vx, *vy, *vz;
    int ns;
    double *np, *mass;
    double cv, bx0, by0, bz0;

    vx = mxGetPr(prhs[0]);
    vy = mxGetPr(prhs[1]);
    vz = mxGetPr(prhs[2]);
    ns = mxGetScalar(prhs[3]);
    np = mxGetPr(prhs[4]);
    mass = mxGetPr(prhs[5]);
    cv = mxGetScalar(prhs[6]);
    bx0 = mxGetScalar(prhs[7]);
    by0 = mxGetScalar(prhs[8]);
    bz0 = mxGetScalar(prhs[9]);
#endif

    /* �o�͍s�������������*/
    plhs[0] = mxCreateDoubleMatrix((mwSize)ns, (mwSize)1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)ns, (mwSize)1, mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *ek, *At;
    ek = mxGetDoubles(plhs[0]);
    At = mxGetDoubles(plhs[1]);
#else
    double *ek, *At;
    ek = mxGetPr(plhs[0]);
    At = mxGetPr(plhs[1]);
#endif
    
    kenergy(ek, At, vx, vy, vz, ns, np, mass, cv, bx0, by0, bz0);
}
