#include "mex.h"
#include<math.h>

#if defined(_OPENMP)
#include<omp.h>
#endif

void calcEJ(double * restrict ejx, double * restrict ejy, double * restrict ejz,
               const double * restrict x, const double * restrict y, 
               const double * restrict vx, const double * restrict vy, const double * restrict vz,
               const int ns, const double * restrict np, const double* restrict q, const double * restrict ex, const double * restrict ey, const double * restrict ez, 
               const int nx, const int ny){

  double xp, yp;
  int ixp, iyp;

  double xw, yw;
  double sf1, sf2, sf3, sf4;
  double p_ex, p_ey, p_ez;
  double p_jx, p_jy, p_jz;

  mxArray *pexavg, *peyavg, *pezavg;

  pexavg = mxCreateDoubleMatrix((mwSize)(nx+2),(mwSize)(ny+2), mxREAL);
  peyavg = mxCreateDoubleMatrix((mwSize)(nx+2),(mwSize)(ny+2), mxREAL);
  pezavg = mxCreateDoubleMatrix((mwSize)(nx+2),(mwSize)(ny+2), mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *exavg, *eyavg, *ezavg;
    exavg = mxGetDoubles(pexavg);
    eyavg = mxGetDoubles(peyavg);
    ezavg = mxGetDoubles(pezavg);
#else
    double *exavg, *eyavg, *ezavg;
    exavg = mxGetPr(pexavg);
    eyavg = mxGetPr(peyavg);
    ezavg = mxGetPr(pezavg);
#endif

 //==============================================================================================
 // Field cancelation to prevent from self-force oscillation
 //==============================================================================================  
  int i_y_x, i_y_xm1, i_ym1_x, i_ym1_xm1;
  
  for(int ix=1;ix<nx+2;ix++){
    for(int iy=1;iy<ny+2;iy++){
        i_y_x     =  iy*(nx+2)    +ix  ; // iy,  ix
        i_y_xm1   =  iy*(nx+2)    +ix-1; // iy,  ix-1
        i_ym1_x   =  (iy-1)*(nx+2)+ix  ; // iy-1,  ix
        i_ym1_xm1 =  (iy-1)*(nx+2)+ix-1; //iy-1, ix-1
        
        exavg[i_y_x] = (                              ex[i_y_xm1] + ex[i_y_x])*0.5;
        eyavg[i_y_x] = (                ey[i_ym1_x]               + ey[i_y_x])*0.5;
        ezavg[i_y_x] = (ez[i_ym1_xm1] + ez[i_ym1_x] + ez[i_y_xm1] + ez[i_y_x])*0.25;
    }
  }
  //--------------------------
  // Main loop
  //--------------------------
  int is, ip;
  int ifirst;
  int ilast;  
  int i1, i2, i3, i4;
  int isidx;

  ifirst=0;
  ilast=0;  

  for(is=0;is<ns;is++){
    ifirst = ilast;
    ilast = ifirst + (int64_T) np[is];
    isidx = is *(nx+2)*(ny+2);

    for(ip=ifirst;ip<ilast;ip++){
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

        p_ex = sf1*exavg[i1] +sf2*exavg[i2] + sf3*exavg[i3] + sf4*exavg[i4];
        p_ey = sf1*eyavg[i1] +sf2*eyavg[i2] + sf3*eyavg[i3] + sf4*eyavg[i4];
        p_ez = sf1*ezavg[i1] +sf2*ezavg[i2] + sf3*ezavg[i3] + sf4*ezavg[i4];

        p_jx = q[is]*vx[ip]*p_ex;
        p_jy = q[is]*vy[ip]*p_ey;
        p_jz = q[is]*vz[ip]*p_ez;

        i1 += isidx;
        i2 += isidx;
        i3 += isidx;
        i4 += isidx;

        ejx[i1] += p_jx*sf1;
        ejx[i2] += p_jx*sf2;
        ejx[i3] += p_jx*sf3;
        ejx[i4] += p_jx*sf4;

        ejy[i1] += p_jy*sf1;
        ejy[i2] += p_jy*sf2;
        ejy[i3] += p_jy*sf3;
        ejy[i4] += p_jy*sf4;

        ejz[i1] += p_jz*sf1;
        ejz[i2] += p_jz*sf2;
        ejz[i3] += p_jz*sf3;
        ejz[i4] += p_jz*sf4;
    }
  }
}

void mexFunction( int nlhs, mxArray *plhs[],       // ??�߂�l??�̌�, ??�߂�l
                  int nrhs, const mxArray *prhs[]) // ??�??�??�??�??�̌�, ??�??�??�??�
{
     /* check for proper number of arguments */
    if(nrhs!=11) {
        mexErrMsgIdAndTxt("c_current:RHS","12 inputs required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("c_current:LHS","3 outputs required.");
    }
    
    /* make sure the first 5 argument are 1d double vector */
    for(int i=0;i<5;i++){
      if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxGetNumberOfDimensions(prhs[i])!=2 ) {
          mexErrMsgIdAndTxt("c_current:particle", "Input multiplier must be a 1d double vector.");
      }
    }

    /* make sure the 6th input argument is scalar */
    if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mxGetNumberOfElements(prhs[5])!=1 ) {
        mexErrMsgIdAndTxt("c_current:ns","6th Input must be a scalar.");
    }

    /* make sure the 7th input argument 1d int vector */
    if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || mxGetNumberOfDimensions(prhs[6])!=2 ) {
        mexErrMsgIdAndTxt("c_current:np", "Input multiplier must be a 1d integer vector.");
    }

    /* make sure the 8th input argument double int vector */
    if( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) || mxGetNumberOfDimensions(prhs[7])!=2 ) {
        mexErrMsgIdAndTxt("c_current:q", "Input multiplier must be a 1d double vector.");
    }

    
    /* make sure the 6th input argument is matrix */
    for(int i=8;i<11;i++){
      if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxGetNumberOfDimensions(prhs[i])!=2) {
          mexErrMsgIdAndTxt("c_current:field size","Input multiplier must be a 2D matrix.");
      }
    }
    /*??�??�??�??�??�̍s??�??�ƃT??�C??�Y??�??�??�擾??�??�??�??�*/    
#if MX_HAS_INTERLEAVED_COMPLEX
    int ns;
    mxDouble *x, *y, *vx, *vy, *vz;
    mxDouble *np, *q;
    mxDouble *ex, *ey, *ez;

    x   = mxGetDoubles(prhs[0]);
    y   = mxGetDoubles(prhs[1]);
    vx = mxGetDoubles(prhs[2]);
    vy = mxGetDoubles(prhs[3]);
    vz = mxGetDoubles(prhs[4]);
    ns = mxGetScalar(prhs[5]);
    np = mxGetDoubles(prhs[6]);
    q = mxGetDoubles(prhs[7]);
    ex = mxGetDoubles(prhs[8]);
    ey = mxGetDoubles(prhs[9]);
    ez = mxGetDoubles(prhs[10]);
#else
    int ns;
    double *x, *y, *vx, *vy, *vz;
    double *np, *q;
    double *ex, *ey, *ez;

    x   = mxGetPr(prhs[0]);
    y   = mxGetPr(prhs[1]);
    vx = mxGetPr(prhs[2]);
    vy = mxGetPr(prhs[3]);
    vz = mxGetPr(prhs[4]);
    ns = mxGetScalar(prhs[5]);
    np = mxGetPr(prhs[6]);
    q = mxGetPr(prhs[7]);
    ex = mxGetPr(prhs[8]);
    ey = mxGetPr(prhs[9]);
    ez = mxGetPr(prhs[10]);
#endif

    /* ??�o??�͍s??�??�??�??�??�??�??�??�??�??�??�??�??�*/
    const mwSize dims[]={mxGetM(prhs[9]),mxGetN(prhs[9]),3};

    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *ejx, *ejy, *ejz;
    ejx = mxGetDoubles(plhs[0]);
    ejy = mxGetDoubles(plhs[1]);
    ejz = mxGetDoubles(plhs[2]);
#else
    double *ejx, *ejy, *ejz;
    ejx = mxGetPr(plhs[0]);
    ejy = mxGetPr(plhs[1]);
    ejz = mxGetPr(plhs[2]);
#endif

    int nx = mxGetM(prhs[9])-2;
    int ny = mxGetN(prhs[9])-2;
    calcEJ(ejx, ejy, ejz,
              x, y, vx, vy, vz,
              ns, np, q,
              ex, ey, ez,
              nx, ny);
}
