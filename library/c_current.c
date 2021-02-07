#include "mex.h"
#include<math.h>

void current(double * restrict jx, double * restrict jy, double * restrict jz, const int nx, const int ny, const double * restrict x, const double * restrict y, const double * restrict vx, const double * restrict vy, const double *restrict vz, const int ns, const double * restrict np, const double* restrict q){

  double xp, yp;
  int ixp, iyp;

  double xw, yw;
  double sf1, sf2, sf3, sf4;
  int i1, i2, i3, i4;
  double p_jx, p_jy, p_jz;
  
  //--------------------------
  // Main loop
  //--------------------------
  int ifirst=0;
  int ilast=0;
  for(int is=0;is<ns;is++){
    ifirst = ilast;
    ilast = ifirst + (int64_T) np[is];

    for(int ip=ifirst;ip<ilast;ip++){
      xp = x[ip] + 1.0;
      yp = y[ip] + 1.0;

      ixp = floor(xp);
      iyp = floor(yp);

      xw = xp - (double) ixp;
      yw = yp - (double) iyp;

      sf3 = xw*yw;
      sf2 = xw-sf3;
      sf4 = yw-sf3;
      sf1 = 1.0-xw-sf4;
      
      i1 = iyp    *(nx+2)+ixp  ; // iyp,  ixp
      i2 = iyp    *(nx+2)+ixp+1; // iyp,  ixp+1
      i3 = (iyp+1)*(nx+2)+ixp+1; // iyp+1,ixp+1
      i4 = (iyp+1)*(nx+2)+ixp  ; // iyp+1,ixp

      p_jx = q[is]*vx[ip];
      p_jy = q[is]*vy[ip];
      p_jz = q[is]*vz[ip];

      jx[i1] += p_jx*sf1;
      jx[i2] += p_jx*sf2;
      jx[i3] += p_jx*sf3;
      jx[i4] += p_jx*sf4;

      jy[i1] += p_jy*sf1;
      jy[i2] += p_jy*sf2;
      jy[i3] += p_jy*sf3;
      jy[i4] += p_jy*sf4;

      jz[i1] += p_jz*sf1;
      jz[i2] += p_jz*sf2;
      jz[i3] += p_jz*sf3;
      jz[i4] += p_jz*sf4;
    }
  }
}

void mexFunction( int nlhs, mxArray *plhs[],       // ?½ß‚ï¿½l?½ÌŒÂï¿½, ?½ß‚ï¿½l
                  int nrhs, const mxArray *prhs[]) // ?½?½?½?½?½ÌŒÂï¿½, ?½?½?½?½
{
     /* check for proper number of arguments */
    if(nrhs!=10) {
        mexErrMsgIdAndTxt("c_current:RHS","8 inputs required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("c_current:LHS","3 outputs required.");
    }
    
    /* make sure the first 5 argument is 1d double vector */
    for(int i=0;i<5;i++){
      if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxGetNumberOfDimensions(prhs[i])!=2 ) {
          mexErrMsgIdAndTxt("c_current:particle", "Input multiplier must be a 1d double vector.");
      }
    }

    /* make sure the 6th input argument is scalar */
    if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mxGetNumberOfElements(prhs[5])!=1 ) {
        mexErrMsgIdAndTxt("c_current:ns","Input multiplier must be a scalar.");
    }

    /* make sure the 7th input argument 1d int vector */
    if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || mxGetNumberOfDimensions(prhs[6])!=2 ) {
        mexErrMsgIdAndTxt("c_current:np", "Input multiplier must be a 1d integer vector.");
    }

    /* make sure the 8th input argument double int vector */
    if( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) || mxGetNumberOfDimensions(prhs[7])!=2 ) {
        mexErrMsgIdAndTxt("c_current:q", "Input multiplier must be a 1d double vector.");
    }

    /* make sure the 6th input argument is scalar */
    for(int i=8;i<10;i++){
      if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxGetNumberOfElements(prhs[i])!=1 ) {
          mexErrMsgIdAndTxt("c_current:field size","Input multiplier must be a scalar.");
      }
    }
    /*?½?½?½?½?½Ìs?½?½ÆƒT?½C?½Y?½?½?½æ“¾?½?½?½?½*/    
#if MX_HAS_INTERLEAVED_COMPLEX
    int ns, nx, ny;
    mxDouble *np; 
    mxDouble *x, *y, *vx, *vy, *vz;
    mxDouble *q;
    
    x = mxGetDoubles(prhs[0]);
    y = mxGetDoubles(prhs[1]);
    vx = mxGetDoubles(prhs[2]);
    vy = mxGetDoubles(prhs[3]);
    vz = mxGetDoubles(prhs[4]);
    ns = mxGetScalar(prhs[5]);
    np = mxGetDoubles(prhs[6]);
    q = mxGetDoubles(prhs[7]);
    nx = mxGetScalar(prhs[8]);
    ny = mxGetScalar(prhs[9]);
#else
    int ns, nx, ny;
    double *np;
    double *x, *y, *vx, *vy, *vz;
    double *q;

    x = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    vx = mxGetPr(prhs[2]);
    vy = mxGetPr(prhs[3]);
    vz = mxGetPr(prhs[4]);
    ns = mxGetScalar(prhs[5]);
    np = mxGetPr(prhs[6]);
    q = mxGetPr(prhs[7]);
    nx = mxGetScalar(prhs[8]);
    ny = mxGetScalar(prhs[9]);
#endif

    /* ?½o?½Ís?½?½?½?½?½?½?½?½?½?½?½?½?½*/
    plhs[0] = mxCreateDoubleMatrix((mwSize)(nx+2), (mwSize)(ny+2), mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)(nx+2), (mwSize)(ny+2), mxREAL);
    plhs[2] = mxCreateDoubleMatrix((mwSize)(nx+2), (mwSize)(ny+2), mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *jx, *jy, *jz;
    jx = mxGetDoubles(plhs[0]);
    jy = mxGetDoubles(plhs[1]);
    jz = mxGetDoubles(plhs[2]);
#else
    double *jx, *jy, *jz;
    jx = mxGetPr(plhs[0]);
    jy = mxGetPr(plhs[1]);
    jz = mxGetPr(plhs[2]);
#endif
    
    current(jx, jy, jz, nx, ny, x, y, vx, vy, vz, ns, np, q);
}
