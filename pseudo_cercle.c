//
//  pseudo_cercle.c
//
//
//  Created by Chantal Oberson Ausoni on 17/09/2015.
//
//

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
    double ux,uy,uz;
    int i,j,k;
    int nber;
    double r,x0,z0;
    
    
    if (argc<4)
    {
        printf("Mesh resolution n m\n");
        exit(EXIT_FAILURE);
    }
    
    if (argc>3)
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
        if(sscanf(argv[3],"%d", &m0) !=1)
        {
            printf("Integer awaited");
            exit(EXIT_FAILURE);
        }
    }
    
    x0 = sin(m0*M_PI/(2*m));
    z0 = cos(m0*M_PI/(2*m))+log(tan(0.5*m0*M_PI/(2*m)));;//point where we want to compute the curvature
    
    fmesh=fopen("pseudo.mesh","w");
    fprintf (fmesh,"MeshVersionFormatted 2\n Dimension\n\t3\n Vertices\n\t%d\n", (m-1)*n+2*m+3);
    //fprintf (fmesh,"MeshVersionFormatted 2\n Dimension\n\t2\n Vertices\n\t%d\n", (m-1)*n);
    
    printf("Make pseudosphere with resolution %d %d\n",n,m);
    printf("Draw circles for point %f %f\n",x0,z0);
    
    for (j=1; j< m; j++)
    {
        for (i=0; i< n; i++)
        {
            x=cos(i*2*M_PI/n)*sin(j*M_PI/(2*m));
            y=sin(i*2*M_PI/n)*sin(j*M_PI/(2*m));
            z=cos(j*M_PI/(2*m))+log(tan(0.5*j*M_PI/(2*m)));
            fprintf(fmesh,"%f\t%f\t%f\t1\n", x, y, z);
        }
    }
    printf("Printed pseudosphere points\n");
    r= x0/sqrt(1-x0*x0);//fabs(tan(z0/x0));
    
    printf("Both circle radii %f %f\n",r,1/r);
    
    for (k=0; k<= m; k++)
    {
        x = r * cos(k*2*M_PI/m);
        y = 0;
        z = r * sin(k*2*M_PI/m)+ z0 + sqrt(r*r-x0*x0);
        ux = 1/r*x0;
        uy = 0;
        uz = 1/r*(-1)*sqrt(r*r-x0*x0);
        xr = ux*ux*(x-x0)+ux*uz*(z-z0)+x0;
        yr = uz*(x-x0)-ux*(z-z0);
        zr = ux*uz*(x-x0)+uz*uz*(z-z0)+z0;
        fprintf(fmesh,"%f\t%f\t%f\t1\n", xr, yr, zr);
        //printf(" Point F: %f\t%f\t%f\t1\n", x, y, z);
        
    }

    printf ("Printed first circle points\n");
    
    
    for (k=0; k<= m; k++)
    {
        x= x0 + 1/(r*r)*x0 + 1/r * cos(k*2*M_PI/m);
        y = 0;
        z= z0 -  1/(r*r)*sqrt(r*r-x0*x0) + 1/r * sin(k*2*M_PI/m);
        fprintf(fmesh,"%f\t%f\t%f\t1\n", x, y, z);
        //printf(" Point F: %f\t%f\t%f\t1\n", x, y, z);
        
    }

    
    printf ("Printed second circle points\n");
    
    fprintf(fmesh,"%f\t%f\t%f\t2\n", x0, 0., z0);
    
    fprintf (fmesh,"Corners\n%d\n%d\n",1,(m-1)*n+2*m+3);
    
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
    
    fprintf (fmesh,"Edges\n%d\n",2*(m+1));
    
    nber = (m-1)*n;
    for (k=1; k<= m; k++)
    {
        fprintf(fmesh,"%d\t%d\t0\n",nber+k,nber+k+1);
    }
    nber = (m-1)*n+(m+1);
    for (k=1; k<= m; k++)
    {
        fprintf(fmesh,"%d\t%d\t0\n",nber+k,nber+k+1);
    }
    fprintf(fmesh,"%d\t%d\t0\n",nber+m+1,nber+1);
    
    printf("Edges printed\n");
    
    fprintf (fmesh,"Ridges\n%d\n",2*(m+1));
    for (k=1; k<= 2*(m+1); k++)
    {
        fprintf(fmesh,"%d\n",k);
    }
    fprintf(fmesh,"End\n");

    printf("Ridges printed\n");
    
    fclose(fmesh);
}
