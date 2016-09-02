//
//  pseudo.c
//  
//
//  Created by Chantal Oberson Ausoni on 16/09/2015.
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
    int n,m;
    double x,y,z;
    int i,j;
    int nber;
    
    
    if (argc<3)
    {
        printf("Mesh resolution n m awaited\n");
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
    fprintf (fmesh,"MeshVersionFormatted 2\n Dimension\n\t3\n Vertices\n\t%d\n", (m-1)*n);
    //fprintf (fmesh,"MeshVersionFormatted 2\n Dimension\n\t2\n Vertices\n\t%d\n", (m-1)*n);
    
    printf("Make pseudosphere with resolution %d %d\n",n,m);
    
    for (j=1; j< m; j++)
    {
        for (i=0; i< n; i++)
        {
            x=cos(i*2*M_PI/n)*sin(j*M_PI/(2*m));
            y=sin(i*2*M_PI/n)*sin(j*M_PI/(2*m));
            z=cos(j*M_PI/(2*m))+log(tan(0.5*j*M_PI/(2*m)));
//            x = (double)i/n;
        //y = (double)j/n;
            fprintf(fmesh,"%f\t%f\t%f\t1\n", x, y, z);
            //fprintf(fmesh,"%f\t%f\t\t1\n", x, y);
            //fprintf(fsol,"%f\n",z);
            
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
    
    fclose(fmesh);
}

