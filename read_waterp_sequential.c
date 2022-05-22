// waterp.c
#include <stdio.h>
#include <math.h>
#include <time.h> // for CPU time
#include <sys/time.h> //for gettimeofday
#include <mpi.h>

#define LENGTH 12000
// global variables
const int maxnum=100000;
double r[maxnum][3][3],rcutsq=1.44,L;
// r(number of molecule, atom 0=O,1=H,2=H, coordinate 0=x,1=y,2=z)

int max(int num1, int num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

int min(int num1, int num2) 
{
    return (num1 > num2 ) ? num2 : num1;
}

double sqr(double a){return a*a;}

double energy12(int i1,int i2){
// ============================
  int m,n,xyz;
  double shift[3],dr[3],mn[3],r6,distsq,dist,ene=0;
  const double sig=0.3166,eps=0.65,eps0=8.85e-12,e=1.602e-19,Na=6.022e23,q[3]={-0.8476,0.4238,0.4238};
  double elst,sig6;
  elst=e*e/(4*3.141593*eps0*1e-9)*Na/1e3,sig6=pow(sig,6);

  // periodic boundary conditions
  for(xyz=0;xyz<=2;xyz++){
   dr[xyz]=r[i1][0][xyz]-r[i2][0][xyz];shift[xyz]=-L*floor(dr[xyz]/L+.5); //round dr[xyz]/L to nearest integer
   dr[xyz]=dr[xyz]+shift[xyz];
  }
  distsq=sqr(dr[0])+sqr(dr[1])+sqr(dr[2]);
  if(distsq<rcutsq){ // calculate energy if within cutoff
    r6=sig6/pow(distsq,3);
    ene=4*eps*r6*(r6-1.); // LJ energy
    for(m=0;m<=2;m++){
      for(n=0;n<=2;n++){
        for(xyz=0;xyz<=2;xyz++) mn[xyz]=r[i1][m][xyz]-r[i2][n][xyz]+shift[xyz];
        dist=sqrt(sqr(mn[0])+sqr(mn[1])+sqr(mn[2]));
        ene=ene+elst*q[m]*q[n]/dist;
    } }
  }
  return ene;
}

main(int argc, char *argv[])
{

  struct timeval start, end;
  clock_t cputime; /* clock_t defined in <time.h> and <sys/types.h> as int */

///////////////////first timer ///////////////////////
cputime = clock();    // assign initial CPU time (IN CPU CLOCKS)
gettimeofday(&start, NULL); // returns structure with time in s and us (microseconds)
 ///////////////////first timer ///////////////////////

  int i,j,natoms,nmol,npairs,me,nproc,k;
  double total,energy=0,Ttime,cpu_total,cpu_sum, max_calc;
  double read_wall_time,read_cpu_time,calc_wall_time,calc_cpu_time,max_read_wall_time,max_cpu_wall_time;
  double before_read_wall_time,before_read_cpu_time,CALC_cpu_time,Bcast_cpu_time,Bcast_wall_time;
  int calc_max, calc_min;

  int minnpair,maxnpair,theo_npair, actual_nparir;
  double mintime, maxtime, imbalance, speedup, overhead;

  FILE *fp;
  char line[LENGTH],nothing[LENGTH],name[20];




 MPI_Init(&argc,&argv); // initialise MPI
 MPI_Comm_size(MPI_COMM_WORLD,&nproc); // return total number of processors
 MPI_Comm_rank(MPI_COMM_WORLD,&me); // return number of this processor, me=0..nproc-1

/////////////////////end first timer ///////////////////////
  cputime= clock()-cputime;      // calculate  cpu clock time as difference of times after-before
  gettimeofday(&end, NULL);
  before_read_wall_time = ((end.tv_sec  - start.tv_sec)+(end.tv_usec - start.tv_usec)/1e6);
  before_read_cpu_time = (float) cputime/CLOCKS_PER_SEC;
/////////////////////end first timer ///////////////////////

/////////////////////begin second timer ///////////////////////
 cputime = clock();    // assign initial CPU time (IN CPU CLOCKS)
 gettimeofday(&start, NULL); // returns structure with time in s and us (microseconds)
////////////////////begin second timer ///////////////////////

for(k=0;k<nproc;k++){
if(me==k){
  fp=fopen("100k.gro", "r"); //opening of file and beginning of reading from HDD
  fgets(line, LENGTH,fp); //skip first line
  fgets(line, LENGTH,fp); sscanf(line,"%i",&natoms);
  nmol= natoms/3; 



  for (i=0;i<nmol;i++){
    for(j=0;j<=2;j++){
      fgets(line, LENGTH,fp);
      sscanf(line, "%s %s %s %lf %lf %lf",nothing,nothing,nothing, &r[i][j][0],&r[i][j][1],&r[i][j][2]);
  } }
  // printf("first line %lf %lf %lf\n",r[0][0][0],r[0][0][1],r[0][0][2]);
  fscanf(fp, "%lf",&L); // read box size


  printf("c: %i Number of molecules %i\n",me,nmol);
  printf("Box size %lf\n",L);
  
   fclose(fp);
}
MPI_Barrier(MPI_COMM_WORLD); // forced waiting of all processes till 

}


/////////////////////end second timer ///////////////////////
  cputime= clock()-cputime;      // calculate  cpu clock time as difference of times after-before
  gettimeofday(&end, NULL);
  read_wall_time = ((end.tv_sec  - start.tv_sec)+(end.tv_usec - start.tv_usec)/1e6);
  read_cpu_time = (float) cputime/CLOCKS_PER_SEC;
/////////////////////end first timer ///////////////////////


////////////////// start fourth timer ////////////
  cputime = clock();    // assign initial CPU time (IN CPU CLOCKS)
  gettimeofday(&start, NULL); // returns structure with time in s and us (microseconds)
////////////////// start fourth timer ////////////

////////////////// calculate energy ////////////////
int low= ( (me* nmol) /nproc) ;
int high=  ((me+1)*nmol)/nproc;

int maxpair,minpair;
long int counter = 0;

  for(i=low;i<high;i++)
  { 

    for(j=i+1;j<nmol;j++)
              {

    energy=energy+energy12(i,j);

    counter = counter + 1;

              }
    

  }
printf("core %i count= %d\n ",me,counter);


////////////////// end fourth timer ////////////
  cputime= clock()-cputime;      // calculate  cpu clock time as difference of times after-before
  gettimeofday(&end, NULL);
  calc_wall_time = ((end.tv_sec  - start.tv_sec)+(end.tv_usec - start.tv_usec)/1e6);
  calc_cpu_time = (float) cputime/CLOCKS_PER_SEC;
////////////////// end fourth timer ////////////



Ttime = calc_wall_time  + read_wall_time;

MPI_Reduce(&calc_wall_time,&max_calc,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD); 
MPI_Allreduce(&energy,&total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 
MPI_Allreduce(&read_cpu_time,&max_cpu_wall_time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 
MPI_Allreduce(&read_wall_time,&max_read_wall_time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 

imbalance = 2.*(maxnpair-minnpair)/(0.+maxnpair+minnpair);


    


  if(me==0){
  printf("Total parallelized energy : %lf \n ",total);
  printf("Energy per molecule %lf \n",total/nmol);

  printf("before read wall time: %lf\n", before_read_wall_time);
  printf("before read CPU  time: %lf\n", before_read_cpu_time);
 printf("max read wall time: %lf\n", max_read_wall_time);
  printf("max read CPU  time: %lf\n", max_cpu_wall_time);
  printf("calc wall time: %lf\n", calc_wall_time);
  printf("calc CPU  time: %lf\n", calc_cpu_time);
  printf("MAX calc wall time: %lf\n", max_calc);
  printf("Total wall  time: %lf\n", Ttime);
  printf("Theoretical time: %lf\n", Ttime/nproc);

  

 }

  //  fclose(fp);


MPI_Finalize(); // finalize MPI peacefully (the system would kill the processes otherwise)



}
