#include "func.h"
double f1(double x) {return 3*exp(x)+sin(x);};
double f2(double x) {return 2*sin(x)+2*cos(x)+exp(x);};
double p1(double x) {return 1;};
double p2(double x) {return 2;};
double q11(double x) {return 1;};
double q12(double x) {return 1;};
double q21(double x) {return 1;};
double q22(double x) {return 3;};

double* func1(double x, double* z)
{
    double* res=(double*)calloc(2,sizeof(double));
    res[0]=z[1];
    res[1]=f1(x)-p1(x)*z[1]-q11(x)*z[0];
    return res;
};

double* func2(double x, double* z)
{
    double* res=(double*)calloc(4,sizeof(double));
    res[0]=z[1];
    res[1]=f1(x)-p1(x)*z[1]-q11(x)*z[0]-q12(x)*z[2];
    res[2]=z[3];
    res[3]=f2(x)-p2(x)*z[3]-q21(x)*z[0]-q22(x)*z[2];
    return res;
};

double func1_acc(double x)
{
    return (1 + sqrt(x) +x)*exp(-x);
}

double* func2_acc(double x)
{
    double* res=(double*)calloc(2,sizeof(double));
    res[0]=exp(x);
    res[1]=sin(x);
    return res;
}