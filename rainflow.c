/*  RAINFLOW $ Revision: 1.1 $ */
/*  by Adam Nieslony, 2009     */

#include <math.h>
#include "mex.h"

/* ++++++++++ BEGIN RF3 [ampl ampl_mean nr_of_cycle] */
/* ++++++++++ Rain flow without time analysis */
void
rf3(mxArray *array_ext, mxArray *hs[]) {
    double *pr, *po, a[16384], ampl, mean;
    int tot_num, index, j, cNr;
    mxArray *array_out;
    
    tot_num = mxGetM(array_ext) * mxGetN(array_ext);
    pr = (double *)mxGetPr(array_ext);
    
    array_out = mxCreateDoubleMatrix(3, tot_num-1, mxREAL);
    po = (double *)mxGetPr(array_out);
    
    j = -1;
    cNr = 1;
    for (index=0; index<tot_num; index++) {
        a[++j]=*pr++;
        while ( (j >= 2) && (fabs(a[j-1]-a[j-2]) <= fabs(a[j]-a[j-1])) ) {
            ampl=fabs( (a[j-1]-a[j-2])/2 );
            switch(j)
{
                case 0: { break; }
                case 1: { break; }
                case 2: {
                    mean=(a[0]+a[1])/2;
                    a[0]=a[1];
                    a[1]=a[2];
                    j=1;
                    if (ampl > 0) {
                        *po++=ampl;
                        *po++=mean;
                        *po++=0.50;
                    }
                    break;
                }
                default: {
                    mean=(a[j-1]+a[j-2])/2;
                    a[j-2]=a[j];
                    j=j-2;
                    if (ampl > 0) {
                        *po++=ampl;
                        *po++=mean;
                        *po++=1.00;
                        cNr++;
                    }
                    break;
                }
            }
        }
    }
    for (index=0; index<j; index++) {
        ampl=fabs(a[index]-a[index+1])/2;
        mean=(a[index]+a[index+1])/2;
        if (ampl > 0){
            *po++=ampl;
            *po++=mean;
            *po++=0.50;
        }
    }
  /* you can free the allocated memeory */
  /* for array_out data                 */
    mxSetN(array_out, tot_num - cNr);
    hs[0]=array_out;
}
/* ++++++++++ END RF3 */

/* ++++++++++ BEGIN RF5 [ampl ampl_mean nr_of_cycle cycle_begin_time cycle_period_time]*/
/* ++++++++++ Rain flow with time analysis */
void
rf5(mxArray *array_ext, mxArray *array_t, mxArray *hs[]) {
    double *pr, *pt, *po, a[16384], t[16384], ampl, mean, period, atime;
    int tot_num, index, j, cNr;
    mxArray *array_out;
    
    tot_num = mxGetM(array_ext) * mxGetN(array_ext);
    pr = (double *)mxGetPr(array_ext);
    pt = (double *)mxGetPr(array_t);
    
    array_out = mxCreateDoubleMatrix(5, tot_num-1, mxREAL);
    po = (double *)mxGetPr(array_out);
    
    j = -1;
    cNr = 1;
    for (index=0; index<tot_num; index++) {
        a[++j]=*pr++;
        t[j]=*pt++;
        while ( (j >= 2) && (fabs(a[j-1]-a[j-2]) <= fabs(a[j]-a[j-1])) ) {
            ampl=fabs( (a[j-1]-a[j-2])/2 );
            switch(j)
{
                case 0: { break; }
                case 1: { break; }
                case 2: {
                    mean=(a[0]+a[1])/2;
                    period=(t[1]-t[0])*2;
                    atime=t[0];
                    a[0]=a[1];
                    a[1]=a[2];
                    t[0]=t[1];
                    t[1]=t[2];
                    j=1;
                    if (ampl > 0) {
                        *po++=ampl;
                        *po++=mean;
                        *po++=0.50;
                        *po++=atime;
                        *po++=period;
                    }
                    break;
                }
                default: {
                    mean=(a[j-1]+a[j-2])/2;
                    period=(t[j-1]-t[j-2])*2;
                    atime=t[j-2];
                    a[j-2]=a[j];
                    t[j-2]=t[j];
                    j=j-2;
                    if (ampl > 0) {
                        *po++=ampl;
                        *po++=mean;
                        *po++=1.00;
                        *po++=atime;
                        *po++=period;
                        cNr++;
                    }
                    break;
                }
            }
        }
    }
    for (index=0; index<j; index++) {
        ampl=fabs(a[index]-a[index+1])/2;
        mean=(a[index]+a[index+1])/2;
        period=(t[index+1]-t[index])*2;
        atime=t[index];
        if (ampl > 0){
            *po++=ampl;
            *po++=mean;
            *po++=0.50;
            *po++=atime;
            *po++=period;
        }
    }
  /* free the memeory !!!*/
    mxSetN(array_out, tot_num - cNr);
    hs[0]=array_out;
}
/* ++++++++++ END RF5 */


/* mexFunction - main function called from MATLAB. */
void
mexFunction( int nlhs,       mxArray *plhs[],
int nrhs, const mxArray *prhs[] )
{
    mxArray *array_in0;
    mxArray *array_in1;
    double *pr, s0, s1, dt;
    int ind;
    
    if (nrhs < 1) {
        mexErrMsgTxt("RAINFLOW requires at least one input argument.");
    } else if (nlhs > 1) {
        mexErrMsgTxt("RAINFLOW requires only one output argument.");
    }
    
    if (mxIsComplex(prhs[0]) || !mxIsDouble(prhs[0])) {
        mexErrMsgTxt("RAINFLOW requires DOUBLE ARRAY as first input argument.");
    } else { array_in0 = (mxArray *)prhs[0]; }
    
    switch(nrhs) {
        case 1: {
            rf3(array_in0, plhs);
            break;
        }
        case 2: {
            if (mxIsComplex(prhs[1]) || !mxIsDouble(prhs[1])) {
                mexErrMsgTxt("RAINFLOW requires two DOUBLE ARRAY input arguments.");
            }
            s0 = mxGetM(prhs[0]) * mxGetN(prhs[0]);
            s1 = mxGetM(prhs[1]) * mxGetN(prhs[1]);
            if (s0 == s1) {
                array_in1 = (mxArray *)prhs[1];
                rf5(array_in0, array_in1, plhs);
            } else if (s1 == 1) {
                pr = (double *)mxGetPr(prhs[1]);
                dt = *pr;
                array_in1 = mxCreateDoubleMatrix(1, s0, mxREAL);
                pr = (double *)mxGetPr(array_in1);
                for (ind=0; ind<s0; ind++) {
                    pr[ind]=ind*dt;
                }
                rf5(array_in0, array_in1, plhs);
                mxDestroyArray(array_in1);
            } else {
                mexErrMsgTxt("RAINFLOW: Time Array size error.");
            }
            break;
        }
        default: {
            mexErrMsgTxt("RAINFLOW: To many input arguments.");
            break;
        }
    }
}