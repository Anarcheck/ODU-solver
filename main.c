#include <stdio.h>
#include "odu.h"
#include "graph.h"

//gcc -Wall -g -o main main.c odu.c func.c -lm

double** calc_y_acc(double a, int points, double h)
{
    double** res=(double**)calloc(points+1,sizeof(double));
    for (int i = 0; i < points+1; i++)
    {
        res[i]=func2_acc(a+i*h);
    }
    return res;
}

int main()
{   
    
    int num_of_eq;
    FILE* fin = fopen("input.txt", "r");
    FILE* fout= fopen("points.txt", "w");
    fscanf(fin,"%d ", &num_of_eq);
    
    if(num_of_eq==1)
    {   
        double a,b,eps,y11,y12;
        int points;
        fscanf(fin,"%lf %lf %lf %d %lf %lf", &a, &b, &eps, &points, &y11, &y12);
        int total_points=0;
        double** answer=solve_boundary_problem1(points, a, b, y11, y12, eps, &total_points);
        double h = (b-a)/points ;
        int step=total_points/points;
        

        for (int i = 0; i < points+1; i++)
        {
            fprintf(fout, "%lf ", a+i*h);
        }
        fprintf(fout, "\n");

        for (int i = 0; i < total_points+1; i+=step)
        {
            fprintf(fout, "%lf ", answer[i][0]);   
        }
        fprintf(fout, "\n");

        for (int i = 0; i < points+1; i++)
        {
            fprintf(fout, "%lf ", func1_acc(a+i*h));
        }
        
        fclose(fout);
        
        for (int i = 0; i < total_points+1; i++)
        {
            free(answer[i]); answer[i]=NULL;
        }
        free(answer); answer=NULL;
        
        if (buildGraph() == 0) {
            printf("Graph successfully build\n");
        }
    }
    
    else if (num_of_eq==2)
    {
        double a,b,eps,y11,y12,y21,y22;
        int points;
        fscanf(fin,"%lf %lf %lf %d %lf %lf %lf %lf", &a, &b, &eps, &points, &y11, &y12, &y21, &y22);
        int total_points=0;
        double** answer=solve_boundary_problem2(points,a,b,y11,y12,y21,y22,eps, &total_points);
        double h = (b-a)/points ;
        int step=total_points/points;
        for (int i = 0; i < points+1; i++)
        {
            fprintf(fout, "%lf ", a+i*h);
        }
        fprintf(fout, "\n");
        
        for (int i = 0; i < total_points+1; i+=step)
        {
            fprintf(fout, "%lf ", answer[i][0]);   
        }
        fprintf(fout, "\n");
        
        for (int i = 0; i < total_points+1; i+=step)
        {
            fprintf(fout, "%lf ", answer[i][2]);   
        }
        fprintf(fout, "\n");
        
        double** y_acc=calc_y_acc(a, points, h);

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < points+1; j++)
            {
                fprintf(fout, "%lf ", y_acc[j][i]);
            }
            fprintf(fout,"\n");
        }
        
        fclose(fout);
        
        for (int i = 0; i < total_points+1; i++)
        {
            free(answer[i]); answer[i]=NULL;
        }
        free(answer); answer=NULL;

        for (int i = 0; i < points+1; i++)
        {
            free(y_acc[i]); y_acc[i]=NULL;
        }
        free(y_acc); y_acc=NULL;
        
        if (buildGraph() == 0) {
            printf("Graph successfully build\n");
        }

    }
    
    else
    {
        printf("Wrong data");
    }
    return 1;               
}