/*The Mandelbrot set is a fractal that is defined as the set of points c
in the complex plane for which the sequence z_{n+1} = z_n^2 + c
with z_0 = 0 does not tend to infinity.*/

/*This code computes an image of the Mandelbrot set.*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define DEBUG 1

#define          X_RESN  1024  /* x resolution */
#define          Y_RESN  1024  /* y resolution */

/* Boundaries of the mandelbrot set */
#define           X_MIN  -2.0
#define           X_MAX   2.0
#define           Y_MIN  -2.0
#define           Y_MAX   2.0

/* More iterations -> more detailed image & higher computational cost */
#define   maxIterations  1000

#define FLOPxI 14
#define MPI_WTIME_IS_GLOBAL 0

typedef struct complextype
{
  float real, imag;
} Compl;

static inline double get_seconds(struct timeval t_ini, struct timeval t_end)
{
  return (t_end.tv_usec - t_ini.tv_usec) / 1E6 +
         (t_end.tv_sec - t_ini.tv_sec);
}

int main(int argc, char *argv[])
{ /* Mandelbrot variables */
  int i, j, k, it = 0, numprocs, rank, Mb, Mbi, Mbim;
  long flopXp;
  Compl   z, c;
  float   lengthsq, temp;
  int *vres, *vresTotal = NULL;

  /* Timestamp variables */
  double tOp0, tOp1, tCo0, tCo1;
  //struct timeval  ti, tf;

  MPI_Init(&argc, &argv);


  tCo0 = MPI_Wtime();

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  tCo1 = MPI_Wtime()-tCo0;


  Mbim = (Y_RESN + numprocs-1)/numprocs; //Redondeo de las filas

  Mbi = rank == numprocs-1? Y_RESN-Mbim*(numprocs-1): Mbim;

  int *res[Mbi], Mbp[numprocs], orden[numprocs];
  long flopXps[numprocs];
  double tOps[numprocs], tCos[numprocs];



  /* Allocate result matrix of Y_RESN x X_RESN */
  vres = (int *) malloc(Mbi * X_RESN * sizeof(int));
  if(!vres)
  { if(rank == 0) { fprintf(stderr, "Error allocating memory\n"); }
    return 1;
  }
  for(i = 0; i < Mbi; i++)
  { res[i] = vres + i*X_RESN; }


  /* Start measuring time */
  tOp0 = MPI_Wtime();

  /* Calculate and draw points */
  for(i = rank*Mbim; i < (rank*Mbim + Mbi); i++)
  { for(j = 0; j < X_RESN; j++)
    { z.real = z.imag = 0.0;
      c.real = X_MIN + j * (X_MAX - X_MIN)/X_RESN;
      c.imag = Y_MAX - i * (Y_MAX - Y_MIN)/Y_RESN;
      k = 0;

      do  /* iterate for pixel color */
      { temp = z.real*z.real - z.imag*z.imag + c.real;  //FLOPxI = 1+1+1+1+1
        z.imag = 2.0*z.real*z.imag + c.imag;            //FLOPxI = 5+1+1+1+1
        z.real = temp;                                  //FLOPxI = 9+1
        lengthsq = z.real*z.real+z.imag*z.imag;         //FLOPxI = 10+1+1+1
        k++; it++;
      } while(lengthsq < 4.0 && k < maxIterations);    //FLOPxI = 13+1
      
      if(k >= maxIterations) res[i%Mbim][j] = 0;
      else { res[i%Mbim][j] = k; }
    }
  } //FLOPxI = 14 

  tOp1 = MPI_Wtime()-tOp0;


  flopXp = FLOPxI * it;

  Mb = Mbi*X_RESN;


  tCo0 = MPI_Wtime();

  //MPI_Gather(&Mb, 1, MPI_INT, Mbp, 1, MPI_INT, 0, MPI_COMM_WORLD);

  for(i=0; i < numprocs; i++)
  {
    Mbp[i] = i == (numprocs-1)? (Y_RESN-Mbi*(numprocs-1))*X_RESN: Mbim*X_RESN;
  }

  tCo1 += MPI_Wtime()-tCo0;


  if(rank == 0)
  { vresTotal = (int *) malloc(Y_RESN * X_RESN * sizeof(int));
    if(!vresTotal)
    { fprintf(stderr, "Error allocating memory\n");
      MPI_Finalize();
      return 1;
    }
    for(i = 0; i < numprocs; i++) { orden[i] = i*Mbim*X_RESN; }
  }


  tCo0 = MPI_Wtime();

  MPI_Gatherv(vres, Mb, MPI_INT, vresTotal, Mbp, orden, MPI_INT, 0, MPI_COMM_WORLD);

  tCo1 += MPI_Wtime()-tCo0;


  //Impresion de resultados:
  MPI_Gather(&tOp1, 1, MPI_DOUBLE, tOps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gather(&tCo1, 1, MPI_DOUBLE, tCos, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gather(&flopXp, 1, MPI_LONG, flopXps, 1, MPI_LONG, 0, MPI_COMM_WORLD);

  if(rank == 0)
  { int *resTotal[Y_RESN];
    long flopTotal = 0, MbpTotal = 0;
    double maxFlops = 0, flops, tOpTotal = 0, tCoTotal = 0, maxtOp = 0, maxtCo = 0;

    for(i = 0; i < Y_RESN; i++)
    { resTotal[i] = vresTotal + i*X_RESN; }

    /* End measuring time */
    //gettimeofday(&tf, NULL);
    for(i = 0; i < numprocs; i++)
    { flops = flopXps[i]/tOps[i];
      tOpTotal += tOps[i];
      tCoTotal += tCos[i];
      flopTotal += flopXps[i];
      MbpTotal += Mbp[i];

      fprintf(stderr, "Proceso nº %d:\n\n", i);
      fprintf(stderr, "· Segundos de computacion:\t%lf\n", tOps[i]);
      fprintf(stderr, "· Segundos de comunicacion:\t%lf\n", tCos[i]);
      fprintf(stderr, "· Segundos en total:\t\t%lf\n", tOps[i]+tCos[i]);
      fprintf(stderr, "· Nº de elementos procesados:\t%d\n", Mbp[i]);
      fprintf(stderr, "· Nº de FLOP hechas:\t\t%ld\n", flopXps[i]);
      fprintf(stderr, "· FLOPS:\t\t\t%lf\n\n", flops);
      
      if(maxFlops < flops) { maxFlops = flops; }
      if(maxtOp < tOps[i]) { maxtOp = tOps[i]; }
      if(maxtCo < tCos[i]) { maxtCo = tCos[i]; }
    }

    flops = flopTotal/maxtOp;

    fprintf(stderr, "\nResumen de todos los procesos:\n\n");
    fprintf(stderr, "· Suma de segundos de computacion:\t\t%lf\n", tOpTotal);
    fprintf(stderr, "· Suma de segundos de comunicacion:\t\t%lf\n", tCoTotal);
    fprintf(stderr, "· Suma de segundos en total:\t\t\t%lf\n", tOpTotal+tCoTotal);
    fprintf(stderr, "· Segundos globales de computacion:\t\t%lf\n", maxtOp);
    fprintf(stderr, "· Segundos globales de comunicacion:\t\t%lf\n", maxtCo);
    fprintf(stderr, "· Segundos globales en total:\t\t\t%lf\n", maxtOp+maxtCo);
    fprintf(stderr, "· Nº de elementos procesados:\t\t\t%ld\n", MbpTotal);
    fprintf(stderr, "· Nº de FLOP hechas:\t\t\t\t%ld\n", flopTotal);
    fprintf(stderr, "· FLOPS:\t\t\t\t\t%lf\n", flops);
    fprintf(stderr, "· Bp = FLOPStotal/(max(FLOPSp)*nprocs) =\t%lf\n", flops/(maxFlops*numprocs));

    /* Print result out */
    if(DEBUG)
    { for(i=0;i<Y_RESN;i++)
      { for(j=0;j<X_RESN;j++)
        { printf("%3d ", resTotal[i][j]); }
        printf("\n");
      }
    }
  }

  free(vres);
  free(vresTotal);

  MPI_Finalize();
  return 0;
}
