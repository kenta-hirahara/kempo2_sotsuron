#include "mex.h"
#include<math.h>

void charge(double * restrict rho, const int nx, const int ny, const double * restrict x, const double * restrict y, const int ns, const double * restrict np, const double* restrict q){

  double xp, yp;
  int ixp, iyp;

  double xw, yw;
  double sf1, sf2, sf3, sf4;
  int i1, i2, i3, i4;

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
      
      rho[i1] += q[is]*sf1;
      rho[i2] += q[is]*sf2;
      rho[i3] += q[is]*sf3;
      rho[i4] += q[is]*sf4;
    }
  }
}

void mexFunction( int nlhs, mxArray *plhs[],       // ?½ß‚ï¿½l?½ÌŒÂï¿½, ?½ß‚ï¿½l
                  int nrhs, const mxArray *prhs[]) // ?½?½?½?½?½ÌŒÂï¿½, ?½?½?½?½
{
     /* check for proper number of arguments */
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("c_current:RHS","8 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("c_current:LHS","3 outputs required.");
    }
    
    /* make sure the first 2 argument is 1d double vector */
    for(int i=0;i<2;i++){
      if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxGetNumberOfDimensions(prhs[i])!=2 ) {
          mexErrMsgIdAndTxt("c_current:particle", "Input multiplier must be a 1d double vector.");
      }
    }

    /* make sure the 3rd input argument is scalar */
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1 ) {
        mexErrMsgIdAndTxt("c_current:ns","Input multiplier must be a scalar.");
    }

    /* make sure the 4th input argument 1d int vector */
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfDimensions(prhs[3])!=2 ) {
        mexErrMsgIdAndTxt("c_current:np", "Input multiplier must be a 1d integer vector.");
    }

    /* make sure the 5th input argument double int vector */
    if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetNumberOfDimensions(prhs[4])!=2 ) {
        mexErrMsgIdAndTxt("c_current:q", "Input multiplier must be a 1d double vector.");
    }

    /* make sure the 6th input argument is scalar */
    for(int i=5;i<7;i++){
      if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxGetNumberOfElements(prhs[i])!=1 ) {
          mexErrMsgIdAndTxt("c_current:field size","Input multiplier must be a scalar.");
      }
    }
    /*?½?½?½?½?½Ìs?½?½ÆƒT?½C?½Y?½?½?½æ“¾?½?½?½?½*/    
#if MX_HAS_INTERLEAVED_COMPLEX
    int ns, nx, ny;
    mxDouble *np;
    mxDouble *x, *y;
    mxDouble *q;
    
    x = mxGetDoubles(prhs[0]);
    y = mxGetDoubles(prhs[1]);
    ns = mxGetScalar(prhs[2]);
    np = mxGetDoubles(prhs[3]);
    q = mxGetDoubles(prhs[4]);
    nx = mxGetScalar(prhs[5]);
    ny = mxGetScalar(prhs[6]);
#else
    int ns, nx, ny;
    double *np;
    double *x, *y;
    double *q;

    x = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    ns = mxGetScalar(prhs[2]);
    np = mxGetPr(prhs[3]);
    q = mxGetPr(prhs[4]);
    nx = mxGetScalar(prhs[5]);
    ny = mxGetScalar(prhs[6]);
#endif

    /* ?½o?½Ís?½?½?½?½?½?½?½?½?½?½?½?½?½*/
    plhs[0] = mxCreateDoubleMatrix((mwSize)(nx+2), (mwSize)(ny+2), mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *rho;
    rho = mxGetDoubles(plhs[0]);
#else
    double *rho;
    rho = mxGetPr(plhs[0]);
#endif
    
    charge(rho, nx, ny, x, y, ns, np, q);
}
