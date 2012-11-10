#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.141592654

void skewcheck(double *x, int *n, int *check)
{
    /* CHECK THE SKEWNESS OF A DISTRIBUTION */
    int i;
    double M=0;
    double sumDist = 0;
    for(i=0; i<*n; i++)
        M += x[i];
    M /= *n;
    for(i=0; i<*n; i++)
        sumDist += pow(x[i]-M,3);
    if(sumDist>0)
        *check = 1;
}

double quantile(double *x, int n, double p)
{
    /* PERCENTILE CALCULATOR (type: 6) */
    int i[2];
    double Q;
    double remain;
    Q = (n+1)*p-1;
    i[0] = floor(Q);
    remain = Q-i[0];
    i[1] = ceil(Q);
    Q = x[i[0]]+(x[i[1]]-x[i[0]])*remain;
    return(Q);
}

double smoothPar(double *x, int n)
{
    /* SMOOTHING PARAMETER */
    int i;
    double M=0;
    double S=0;
    double Q[2];
    for(i=0; i<n; i++)
        M += x[i];
    M /= n;
    for(i=0; i<n; i++)
        S += pow(x[i]-M,2);
    S /= (n-1);
    S = pow(S,0.5);
    Q[0] = quantile(x,n,0.25);
    Q[1] = quantile(x,n,0.75);
    Q[0] = (Q[1]-Q[0])/1.349;
    if(S < Q[0])
        Q[0] = S;
    Q[0] = (0.51/pow(n,0.2))*Q[0];
    return(Q[0]);
}

double kernel(double *x, int n, double h, int j)
{
    /* GAUSSIAN KERNEL */
    int i;
    double k=0;
    double z[n];
    for(i=0; i<n; i++) {
        z[i] = (x[i]-x[j])/h;
        k  += (1/pow(2*PI,0.5))*exp(-pow(z[i],2)/2);
    }
    k = (1/(n*h))*k;
    return(k);
}

void kernestim(double *x, int *n, double *smoothing, double *value)
{
    /* KERNEL DENSITY ESTIMATOR */
    int i;
    int maxD=0;
    double D[2];
    if(*smoothing==0)
        *smoothing = smoothPar(x,*n);
    D[0] = kernel(x,*n,*smoothing,0);
    for(i=1; i<*n; i++) {
        D[1] = kernel(x,*n,*smoothing,i);
        if(D[1]>D[0]) {
            D[0] = D[1];
            maxD = i;
        }
    }
    *value = x[maxD];
}

void histestim(double *x, int *n, double *smoothing, double *value)
{
    /* HISTOGRAM DENSITY ESTIMATOR */
    int i,j;
    int count;
    int peakCount=0;
    double adding=0;
    double lim[2];
    double RSS[2];
    double m;
    RSS[0] = 0;
    RSS[1] = 0;
    if(*smoothing==0) {
        *smoothing = quantile(x,*n,0.975);
        *smoothing -= quantile(x,*n,0.025);
        *smoothing /= pow(*n,0.5);
    }
    *smoothing *= 0.5;
    for(i=0; i<*n; i++) {
        lim[0] = x[i]-(*smoothing);
        lim[1] = x[i]+(*smoothing);
        adding = 0;
        count = 0;
        for(j=0; j<*n; j++) {
            if((x[j]>=lim[0]) & (x[j]<=lim[1])) {
                count++;
                adding += x[j];
            }
        }
        if(count > peakCount) {
            *value = adding/count;
            peakCount = count;
        } else {
            if(count == peakCount) {
                m = adding/count;
                for(j=0; j<*n; j++) {
                    RSS[0] += pow(x[j]-(*value),2);
                    RSS[1] += pow(x[j]-m,2);
                }
                if(RSS[1] < RSS[0]) {
                    *value = m;
                    peakCount = count;
                    RSS[0] = 0;
                    RSS[1] = 0;
                }
            }
        }
    }
}

void distestim(double *x, int *n, double *smoothing, double *value)
{
    /* DISTANCE DENSITY ESTIMATOR
    Around each data point (pivot), an interval [x-h/2,x+h/2] is builded, where h is the
    smoothing parameter. The function searches the interval (bin) with the higher data frequency.
    The output value is the weighted average of the values into the selected bin, in which each
    observation is weighted on the basis of the distance from the pivot. If bins with equal
    densities are found, the bin presenting the smallest deviance from the pivot is chosen.
    */
    int i,j;            // Iterative indices
    int count;          // Frequency of the current bin
    int peakCount=0;    // Frequency of the selected bin
    double lim[2];      // Lower and upper bounds of the current bin
    double range;       // Data range within the current bin
    double peakRange=0; // Data range within the selected bin
    double M;           // Weighted average of the current bin
    double D[*n];       // Distance between data and M into the selected bin
    double sumD=0;      // Sum of the distances D
    double xBin[*n];    // Data of the current bin
    if(*smoothing==0) {
        *smoothing = quantile(x,*n,0.975);
        *smoothing -= quantile(x,*n,0.025);
        *smoothing /= pow(*n,0.5);
    }
    *smoothing *= 0.5;
    for(i=0; i<*n; i++) {
        lim[0] = x[i]-(*smoothing);
        lim[1] = x[i]+(*smoothing);
        sumD = 0;
        count = 0;
        for(j=0; j<*n; j++) {
            if((x[j]>=lim[0]) & (x[j]<=lim[1])) {
                xBin[count] = x[j];
                D[count] = pow(x[j]-x[i],2);
                if(D[count] < 0)
                    D[count] *= -1;
                sumD += D[count];
                count++;
            }
        }
        if(count >= peakCount) {
            range = sumD;
            for(j=0; j<count; j++)
                D[j] = 1-D[j]/sumD;
            sumD = 0;
            M = 0;
            for(j=0; j<count; j++) {
                M += xBin[j]*D[j];
                sumD += D[j];
            }
            M /= sumD;
            if((count > peakCount) || (range < peakRange)){
                *value = M;
                peakCount = count;
                peakRange = range;
            }
        }
    }
}
