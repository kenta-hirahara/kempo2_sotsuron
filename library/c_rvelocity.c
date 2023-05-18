#include "mex.h"
#include<math.h>

#if defined(_OPENMP)
#include<omp.h>
#endif

void rvelocity(double * restrict vxo, double * restrict vyo, double * restrict vzo,
               const double * restrict x, const double * restrict y, 
               const double * restrict vxi, const double * restrict vyi, const double * restrict vzi,
               const int ns, const double * restrict np, const double* restrict qm, const double cv,
               const double * restrict ex, const double * restrict ey, const double * restrict ez, 
               const double * restrict bx, const double * restrict by, const double * restrict bz, 
               const int nx, const int ny){

  double xp, yp;
  int ixp, iyp;

  double xw, yw;
  double sf1, sf2, sf3, sf4;
  double p_ex, p_ey, p_ez, p_bx, p_by, p_bz;

  double boris, gamma;
  double ux, uy, uz, ux_r, uy_r, uz_r, ux_f, uy_f, uz_f;
  const double csq = cv*cv;

  mxArray *pexavg, *peyavg, *pezavg;
  mxArray *pbxavg, *pbyavg, *pbzavg;

  pexavg = mxCreateDoubleMatrix((mwSize)(nx+2),(mwSize)(ny+2), mxREAL);
  peyavg = mxCreateDoubleMatrix((mwSize)(nx+2),(mwSize)(ny+2), mxREAL);
  pezavg = mxCreateDoubleMatrix((mwSize)(nx+2),(mwSize)(ny+2), mxREAL);
  pbxavg = mxCreateDoubleMatrix((mwSize)(nx+2),(mwSize)(ny+2), mxREAL);
  pbyavg = mxCreateDoubleMatrix((mwSize)(nx+2),(mwSize)(ny+2), mxREAL);
  pbzavg = mxCreateDoubleMatrix((mwSize)(nx+2),(mwSize)(ny+2), mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *exavg, *eyavg, *ezavg;
    mxDouble *bxavg, *byavg, *bzavg;
    exavg = mxGetDoubles(pexavg);
    eyavg = mxGetDoubles(peyavg);
    ezavg = mxGetDoubles(pezavg);
    bxavg = mxGetDoubles(pbxavg);
    byavg = mxGetDoubles(pbyavg);
    bzavg = mxGetDoubles(pbzavg);
#else
    double *exavg, *eyavg, *ezavg;
    double *bxavg, *byavg, *bzavg;
    exavg = mxGetPr(pexavg);
    eyavg = mxGetPr(peyavg);
    ezavg = mxGetPr(pezavg);
    bxavg = mxGetPr(pbxavg);
    byavg = mxGetPr(pbyavg);
    bzavg = mxGetPr(pbzavg);
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
        bxavg[i_y_x] = (                              bx[i_y_xm1] + bx[i_y_x])*0.5;
        byavg[i_y_x] = (                by[i_ym1_x]               + by[i_y_x])*0.5;
        bzavg[i_y_x] = (bz[i_ym1_xm1] + bz[i_ym1_x] + bz[i_y_xm1] + bz[i_y_x])*0.25;
    }
  }
  //--------------------------
  // Main loop
  //--------------------------
  int is, ip;
  int ifirst;
  int ilast;  
  int i1, i2, i3, i4;

  #pragma omp parallel private(is, ifirst, ilast)
  {
    ifirst=0;
    ilast=0;  
  
    for(is=0;is<ns;is++){
      ifirst = ilast;
      ilast = ifirst + (int64_T) np[is];

      #pragma omp for simd nowait safelen(8) linear(ip:1) \
      private(xp,yp,ixp,iyp,xw,yw),\
      private(sf1,sf2,sf3,sf4,i1,i2,i3,i4), \
      private(p_ex,p_ey,p_ez,p_bx,p_by,p_bz), \
      private(gamma,boris), \
      private(ux,uy,uz,ux_r,uy_r,uz_r,ux_f,uy_f,uz_f)
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

        sf1*=qm[is];
        sf2*=qm[is];
        sf3*=qm[is];
        sf4*=qm[is];
        
        i1 = iyp    *(nx+2)+ixp  ; // iyp,  ixp
        i2 = iyp    *(nx+2)+ixp+1; // iyp,  ixp+1
        i3 = (iyp+1)*(nx+2)+ixp+1; // iyp+1,ixp+1
        i4 = (iyp+1)*(nx+2)+ixp  ; // iyp+1,ixp

        p_ex = sf1*exavg[i1] +sf2*exavg[i2] + sf3*exavg[i3] + sf4*exavg[i4];
        p_ey = sf1*eyavg[i1] +sf2*eyavg[i2] + sf3*eyavg[i3] + sf4*eyavg[i4];
        p_ez = sf1*ezavg[i1] +sf2*ezavg[i2] + sf3*ezavg[i3] + sf4*ezavg[i4];
        p_bx = sf1*bxavg[i1] +sf2*bxavg[i2] + sf3*bxavg[i3] + sf4*bxavg[i4];
        p_by = sf1*byavg[i1] +sf2*byavg[i2] + sf3*byavg[i3] + sf4*byavg[i4];
        p_bz = sf1*bzavg[i1] +sf2*bzavg[i2] + sf3*bzavg[i3] + sf4*bzavg[i4];

        // STEP1: derive Lorlentz factor
        gamma = cv/sqrt(csq - vxi[ip]*vxi[ip] - vyi[ip]*vyi[ip] - vzi[ip]*vzi[ip]);
        // STEP2: Obtain momentum with accerelation by electric field for half time step
        ux = vxi[ip]*gamma + p_ex;
        uy = vyi[ip]*gamma + p_ey;
        uz = vzi[ip]*gamma + p_ez;

        gamma = cv/sqrt(csq + ux*ux + uy*uy + uz*uz);
        // STEP4:B-field correction by relativistic factor
        p_bx *= gamma;
        p_by *= gamma;
        p_bz *= gamma;

        boris = 2.0/(1.0 + p_bx*p_bx + p_by*p_by + p_bz*p_bz);

        // STEP5:Halfstep rotation by B-field
        ux_r = ux + uy*p_bz - uz*p_by;
        uy_r = uy + uz*p_bx - ux*p_bz;
        uz_r = uz + ux*p_by - uy*p_bx;

        // STEP6:fullstep rotation correction by Boris factor and E-field accerelation
        ux_f = ux + boris*(uy_r*p_bz - uz_r*p_by) + p_ex;
        uy_f = uy + boris*(uz_r*p_bx - ux_r*p_bz) + p_ey;
        uz_f = uz + boris*(ux_r*p_by - uy_r*p_bx) + p_ez;

        gamma = cv/sqrt(csq + ux_f*ux_f + uy_f*uy_f + uz_f*uz_f);

        // STEP7:Momentum to velocity
        vxo[ip] = ux_f*gamma;
        vyo[ip] = uy_f*gamma;
        vzo[ip] = uz_f*gamma;
      }
    }
  }
}

void mexFunction( int nlhs, mxArray *plhs[],       // ?½ß‚ï¿½l?½ÌŒÂï¿½, ?½ß‚ï¿½l
                  int nrhs, const mxArray *prhs[]) // ?½?½?½?½?½ÌŒÂï¿½, ?½?½?½?½
{
     /* check for proper number of arguments */
    if(nrhs!=15) {
        mexErrMsgIdAndTxt("c_current:RHS","8 inputs required.");
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

    /* make sure the 9th input argument is scalar */
    if( !mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) || mxGetNumberOfElements(prhs[8])!=1 ) {
        mexErrMsgIdAndTxt("c_current:cv","9th Input cv must be a scalar.");
    }
    
    /* make sure the 6th input argument is matrix */
    for(int i=9;i<15;i++){
      if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxGetNumberOfDimensions(prhs[i])!=2) {
          mexErrMsgIdAndTxt("c_current:field size","Input multiplier must be a 2D matrix.");
      }
    }
    /*?½?½?½?½?½Ìs?½?½ÆƒT?½C?½Y?½?½?½æ“¾?½?½?½?½*/    
#if MX_HAS_INTERLEAVED_COMPLEX
    int ns;
    mxDouble *x, *y, *vxi, *vyi, *vzi;
    mxDouble *np, *qm;
    double cv;
    mxDouble *ex, *ey, *ez;
    mxDouble *bx, *by, *bz;

    x   = mxGetDoubles(prhs[0]);
    y   = mxGetDoubles(prhs[1]);
    vxi = mxGetDoubles(prhs[2]);
    vyi = mxGetDoubles(prhs[3]);
    vzi = mxGetDoubles(prhs[4]);
    ns = mxGetScalar(prhs[5]);
    np = mxGetDoubles(prhs[6]);
    qm = mxGetDoubles(prhs[7]);
    cv = mxGetScalar(prhs[8]);
    ex = mxGetDoubles(prhs[9]);
    ey = mxGetDoubles(prhs[10]);
    ez = mxGetDoubles(prhs[11]);
    bx = mxGetDoubles(prhs[12]);
    by = mxGetDoubles(prhs[13]);
    bz = mxGetDoubles(prhs[14]);
#else
    int ns;
    double *x, *y, *vxi, *vyi, *vzi;
    double *np, *qm;
    double cv;
    double *ex, *ey, *ez;
    double *bx, *by, *bz;

    x   = mxGetPr(prhs[0]);
    y   = mxGetPr(prhs[1]);
    vxi = mxGetPr(prhs[2]);
    vyi = mxGetPr(prhs[3]);
    vzi = mxGetPr(prhs[4]);
    ns = mxGetScalar(prhs[5]);
    np = mxGetPr(prhs[6]);
    qm = mxGetPr(prhs[7]);
    cv = mxGetScalar(prhs[8]);
    ex = mxGetPr(prhs[9]);
    ey = mxGetPr(prhs[10]);
    ez = mxGetPr(prhs[11]);
    bx = mxGetPr(prhs[12]);
    by = mxGetPr(prhs[13]);
    bz = mxGetPr(prhs[14]);
#endif

    /* ?½o?½Ís?½?½?½?½?½?½?½?½?½?½?½?½?½*/
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[2]), mxGetN(prhs[2]), mxREAL);
    plhs[1] = mxCreateDoubleMatrix(mxGetM(prhs[3]), mxGetN(prhs[3]), mxREAL);
    plhs[2] = mxCreateDoubleMatrix(mxGetM(prhs[4]), mxGetN(prhs[4]), mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *vxo, *vyo, *vzo;
    vxo = mxGetDoubles(plhs[0]);
    vyo = mxGetDoubles(plhs[1]);
    vzo = mxGetDoubles(plhs[2]);
#else
    double *vxo, *vyo, *vzo;
    vxo = mxGetPr(plhs[0]);
    vyo = mxGetPr(plhs[1]);
    vzo = mxGetPr(plhs[2]);
#endif

    int nx = mxGetM(prhs[9])-2;
    int ny = mxGetN(prhs[9])-2;
    rvelocity(vxo, vyo, vzo,
              x, y, vxi, vyi, vzi,
              ns, np, qm, cv,
              ex, ey, ez,
              bx, by, bz,
              nx, ny);
}
