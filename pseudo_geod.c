//
//  pseudo_geod.c
//
//
//  Created by Chantal Oberson Ausoni on 17/09/2015.
//
// generating pseudosphere mesh, with geodesic in a point
// parameters n and m are the resolution of the mesh, m0 determines where the point lies
// output pseudo.mesh
// calling function with: ./a.out 20 20 13
// in medit: g to show principal circles

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[])
{
    FILE* fmesh;
    //FILE* fsol;
    int n,m,m0;
    double x,y,z;
    double xr,yr,zr;
    double p;
    int i,j,k;
    int nber;
    double r,x0,z0;
    
    
    if (argc<3)
    {
        printf("Mesh resolution n m awaited and row of the reference point m0\n");
        exit(EXIT_FAILURE);
    }
    
    if (argc>2)
    {
        if(sscanf(argv[1],"%d", &n) !=1)
        {
            printf("Integer awaited");
            exit(EXIT_FAILURE);
        }
        if(sscanf(argv[2],"%d", &m) !=1)
        {
            printf("Integer awaited");
            exit(EXIT_FAILURE);
        }
    }
    
    
    
    fmesh=fopen("pseudo.mesh","w");
    fprintf (fmesh,"MeshVersionFormatted 2\n Dimension\n\t3\n Vertices\n\t%d\n", (m-1)*n+m);
    //fprintf (fmesh,"MeshVersionFormatted 2\n Dimension\n\t2\n Vertices\n\t%d\n", (m-1)*n);
    
    printf("Make pseudosphere with resolution %d %d\n",n,m);
    
    for (j=1; j< m; j++)
    {
        for (i=0; i< n; i++)
        {
            r=j/(double)m;
            //r= 1/(double)(j*m);
            p=sqrt(1-r*r);
            x=r*cos(i*2*M_PI/n);
            z=r*sin(i*2*M_PI/n);
            y=log(1/r*(1+p))-p;
            fprintf(fmesh,"%f\t%f\t%f\t1\n", x, y, z);
        }
    }
    printf("Printed pseudosphere points\n");
    
    for (k=1; k< m; k++)
    {
        x = 5*M_PI * cos(k*M_PI/m)- 5*M_PI;
        //y = 1 + 3*M_PI * sin(k*M_PI/m);//essai avec 1+ ou test avec if?
        y = 5*M_PI * sin(k*M_PI/m);
        if (y>1)
        {
        r= 1/y; //1/y;   //CHANGE or not?
        p=sqrt(1-r*r);
        xr = r*cos(x);
        zr = r*sin(x);
        yr = log(1/r*(1+p))-p;
        fprintf(fmesh,"%f\t%f\t%f\t1\n", xr, yr, zr);
        printf(" P : %f %f\t%f\t%f\t%f\n",p, x, y,r,yr);
        }
        
    }


    
    fprintf (fmesh,"Triangles\n%d\n",2*n*(m-2));
    for (j=1; j< m-1; j++)
    {
        for (i=1; i< n; i++)
        {
            nber = i+(j-1)*n;
            fprintf(fmesh,"%d\t%d\t%d\t1\n",nber,nber+n,nber+n+1);
            fprintf(fmesh,"%d\t%d\t%d\t1\n",nber,nber+1,nber+1+n);
        }
        
        nber = 1+(j-1)*n;
        fprintf(fmesh,"%d\t%d\t%d\t1\n",nber,nber+2*n-1,nber+n);
        fprintf(fmesh,"%d\t%d\t%d\t1\n",nber,nber+n-1,nber+2*n-1);
        
    }
    
    printf("Triangles printed \n");
    
    fprintf (fmesh,"Edges\n%d\n",m-1);
    
    nber = (m-1)*n;
    for (k=1; k<m; k++)
    {
        fprintf(fmesh,"%d\t%d\t0\n",nber+k,nber+k+1);
    }
    
    fprintf (fmesh,"Ridges\n%d\n",m-1);
    for (k=1; k<m; k++)
    {
        fprintf(fmesh,"%d\n",k);
    }

    
    fclose(fmesh);
}
