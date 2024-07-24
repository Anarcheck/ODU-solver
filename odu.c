#include "odu.h"
#include <stdio.h>

#define EPS 1E-5
double* calcK(double x, double* z, double h, int s)
{
   int n=2*s;
   double *res;
   if (s==1)
   {
      res=func1(x,z);
      for (int i = 0; i < n; i++)
      {
         res[i]*=h;
      }
      
   }
   else if (s==2)
   {
      res=func2(x,z);
      for (int i = 0; i < n; i++)
      {
         res[i]*=h;
      }
   }
   return res;
}

double* calculateNextZ(double x, double* z, double h, int s)
{
      int n=2*s;
      double* new_z = (double*)calloc(n, sizeof(double));
      double* k1=calcK(x,z,h,s);
      
      for (int i = 0; i < n; i++) new_z[i] = z[i] + 0.5*k1[i]; 
      
      double* k2=calcK(x +0.5*h, new_z, h, s);
      
      for (int i = 0; i < n; i++) new_z[i] = z[i]+0.5*k2[i]; 
      
      double* k3=calcK(x+0.5*h, new_z, h, s);
      
      for (int i = 0; i < n; i++) new_z[i] = z[i]+k3[i]; 
      
      double* k4=calcK(x+h, new_z, h ,s);
      for (int i = 0; i < n; i++)
      {
         new_z[i]=z[i] + (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6 ;
      }

      free(k1); k1=NULL;
      free(k2); k2=NULL;
      free(k3); k3=NULL;
      free(k4); k4=NULL;
      return new_z;

}

double** solveRungeKutta(int points, double x, double* z0, double h, int s)
{
   double** res=(double**)calloc(points+1, sizeof(double));
   res[0]=z0;
   x+=h;
   for(int i=1; i<points+1; i++)
   {
      res[i]=calculateNextZ(x, res[i-1], h, s);
      x+=h;
   }
   return res;
}

double GaussSystem(double b1, double b2, double b3, double h1)
{
   double h2 = h1/2;
   double h3 = h1/4;
   double matrix[3][3] = 
   {
        {exp(h1), exp(2 * h1), exp(3 * h1)},
        {exp(h2), exp(2 * h2), exp(3 * h2)},
        {exp(h3), exp(2 * h3), exp(3 * h3)}
   };
   double coef1 = matrix[1][0] / matrix[0][0];
   double coef2 = matrix[2][0] / matrix[0][0];
   matrix[1][0] = matrix[1][0] - coef1 * matrix[0][0];
   matrix[2][0] = matrix[2][0] - coef2 * matrix[0][0];
   matrix[1][1] = matrix[1][1] - coef1 * matrix[0][1];
   matrix[2][1] = matrix[2][1] - coef2 * matrix[0][1];
   matrix[1][2] = matrix[1][2] - coef1 * matrix[0][2];
   matrix[2][2] = matrix[2][2] - coef2 * matrix[0][2];
   b2 = b2 - coef1 * b1;
   b3 = b3 - coef2 * b1;
   coef1 = matrix[2][1] / matrix[1][1];
   matrix[2][1] = matrix[2][1] - coef1 * matrix[1][1];
   matrix[2][2] = matrix[2][2] - coef1 * matrix[1][2];
   b3 = b3 - coef1 * b2;
   double C = b3 / matrix[2][2];
   double B = (b2 - matrix[1][2] * C) / matrix[1][1];
   double A = (b1 - matrix[0][2] * C - matrix[0][1] * B) / matrix[0][0];
   return A + B + C;
}

double** make_approximation(double** z1, double** z2, double** z3, int points, double h, int s)
{
   int n=2*s;
   double** res=(double**)calloc(points+1, sizeof(double));
   for (int i = 0; i < points+1; i++)
   {
      res[i]=(double*)calloc(n, sizeof(double));
      for (int j = 0; j < n; j++)
      {
         res[i][j]=GaussSystem(z1[i][j], z2[2*i][j], z3[4*i][j], h);
         
      }
   }
   return res;
}

double calc_eps(double** a1, double** a2, int points)
{
   double s;
   for (int i = 0; i < points+1; i++)
   {
      s+=pow(a1[i][0]-a2[2*i][0], 2);
   }
   double res=sqrt(s)/points;
   return  res;
   
}

double** approximate1(int points, double x, double y11, double t, double h, double eps, int* total_points)
{
   double epsilon=1;
   double z0[2]={y11, t};
   double** z1 = solveRungeKutta(points, x, z0, h, 1);
   double** z2 = solveRungeKutta(2*points, x, z0, h/2, 1);
   double** z3 = solveRungeKutta(4*points, x, z0, h/4, 1);
   double** z4 = solveRungeKutta(8*points, x, z0, h/8, 1);
   double** approx1 = make_approximation(z1, z2, z3, points, h, 1);
   double** approx2 = make_approximation(z2, z3, z4,  2 * points, h/2, 1);

   do
   {
      epsilon=calc_eps(approx1,approx2,points);
      approx1=approx2;
      z2=z3;
      z3=z4;
      points*=2;
      h/=2;
      z4 = solveRungeKutta(8*points, x, z0, h/8, 1);
   } while (epsilon>eps);
   
   *total_points=points;
   return approx1;
   /* while(1)
   {
      epsilon = calc_eps(approx1,approx2,points);
      if(epsilon<eps)
      {
         *total_points=2*points;
         break;
      }
      approx1=approx2;
      z2=z3;
      z3=z4;
      points*=2;
      h/=2;
      z4 = solveRungeKutta(8*points, x, z0, h/8, 1);
      approx2 = make_approximation(z2, z3, z4, 2*points, h/2, 1);
   }
   return approx2; */
}
double** approximate2(int points, double x, double y11, double y21, double t1, double t2, double h, double eps, int* total_points)
{
   double epsilon=100;
   double z0[4]={y11, t1, y21, t2};
   double** z1 = solveRungeKutta(points, x, z0, h, 2);
   double** z2 = solveRungeKutta(2*points, x, z0, h/2, 2);
   double** z3 = solveRungeKutta(4*points, x, z0, h/4, 2);
   double** z4 = solveRungeKutta(8*points, x, z0, h/8, 2);
   double** approx1 = make_approximation(z1, z2, z3, points, h, 2);
   double** approx2 = make_approximation(z2, z3, z4,  2 * points, h/2, 2);

   do
   {
      epsilon=calc_eps(approx1,approx2,points);
      approx1=approx2;
      z2=z3;
      z3=z4;
      points*=2;
      h/=2;
      z4 = solveRungeKutta(8*points, x, z0, h/8, 2);
   } while (epsilon>eps);
   
   *total_points=points;
   return approx1;
}

double** solve_boundary_problem1(int points, double a, double b, double y11, double y12, double eps, int* total_points)
{
   double t_left, t_right, t_mid;
   t_left=-50;
   t_right=50;
   double h = (b-a)/points;
   int total_points1, total_points2;
   
   while(t_right-t_left > EPS)
   {
      double** res_l=approximate1(points, a, y11, t_left, h, eps, &total_points1);
      t_mid=(t_right+t_left)/2;
      double** res_m=approximate1(points, a, y11, t_mid, h, eps, &total_points2);
      if (fabs(res_m[total_points2][0]-y12)<EPS)
      {
         *total_points=total_points2;
         return res_m;
      }
      else
      {
         if ( (res_m[total_points2][0]-y12)*(res_l[total_points1][0]-y12)>0)
         {
            t_left=t_mid;
         }
         else t_right=t_mid;
      }
   }
   double** res_m=approximate1(points, a, y11, t_mid, h, eps, total_points);
  
   return res_m;
}

double** solve_boundary_problem2(int points, double a, double b, double y11, double y12, double y21, double y22, double eps, int* total_points)
{
   double t1=1;
   double t2=1;
   double h = (b-a)/points;
   double** res=approximate2(points, a, y11, y21, t1, t2, h, eps, total_points);
   return res;
}