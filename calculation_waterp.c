// waterp.c
#include <stdio.h>
#include <math.h>
#include <time.h> // for CPU time
#include <sys/time.h> //for gettimeofday
#include <mpi.h>

#define LENGTH 1000
// global variables
const long int maxnum=1000000;
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

  int i,j,natoms,nmol,npairs,me,nproc;
  double total,energy=0,Ttime,cpu_total,cpu_sum;
  double read_wall_time,read_cpu_time,calc_wall_time,calc_cpu_time;
  double before_read_wall_time,before_read_cpu_time,CALC_cpu_time,Bcast_cpu_time,Bcast_wall_time;
  int calc_max, calc_min;

  int minnpair,maxnpair,theo_npair, actual_nparir;
  double mintime, maxtime, imbalance, speedup, overhead,max_calc;

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
if(me==0){
  printf("Program to calculate energy of water\n");
  printf("Input NAME of configuration file\n");
  scanf("%s",name); // reading of filename from keyboard 

  fp=fopen(name, "r"); //opening of file and beginning of reading from HDD
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


  printf("Number of molecules %i\n",nmol);
  printf("Box size %lf\n",L);
  
}


/////////////////////end second timer ///////////////////////
  cputime= clock()-cputime;      // calculate  cpu clock time as difference of times after-before
  gettimeofday(&end, NULL);
  read_wall_time = ((end.tv_sec  - start.tv_sec)+(end.tv_usec - start.tv_usec)/1e6);
  read_cpu_time = (float) cputime/CLOCKS_PER_SEC;
/////////////////////end first timer ///////////////////////


/////////////////////begin third timer ///////////////////////
 cputime = clock();    // assign initial CPU time (IN CPU CLOCKS)
 gettimeofday(&start, NULL); // returns structure with time in s and us (microseconds)
////////////////////begin third timer ///////////////////////
MPI_Bcast(&nmol,100,MPI_INT,0,MPI_COMM_WORLD); //broadcast name of file from thread 0 to all threads
MPI_Bcast(&r,maxnum*3*3,MPI_DOUBLE,0,MPI_COMM_WORLD); //broadcast name of file from thread 0 to all threads
MPI_Bcast(&L,100,MPI_DOUBLE,0,MPI_COMM_WORLD); //broadcast name of file from thread 0 to all threads

/////////////////////end third timer ///////////////////////
  cputime= clock()-cputime;      // calculate  cpu clock time as difference of times after-before
  gettimeofday(&end, NULL);
  Bcast_wall_time = ((end.tv_sec  - start.tv_sec)+(end.tv_usec - start.tv_usec)/1e6);
  Bcast_cpu_time = (float) cputime/CLOCKS_PER_SEC;
/////////////////////end third timer ///////////////////////


////////////////// start fourth timer ////////////
  cputime = clock();    // assign initial CPU time (IN CPU CLOCKS)
  gettimeofday(&start, NULL); // returns structure with time in s and us (microseconds)
////////////////// start fourth timer ////////////

////////////////// calculate energy ////////////////
int low= ( (me* nmol) /nproc) ;
int high=  ((me+1)*nmol)/nproc;

int maxpair,minpair;
long int counter = 0;

  for(i=me;i<nmol;i=i+nproc)
  { 

    for(j=i+1;j<nmol;j++)
              {

    energy=energy+energy12(i,j);

    counter = counter + 1;
    // if (count == (nmol*.25)) { me=2;}
    // else if (count == (nmol*.5)) { me=1;}

    // else if (count == (nmol*.75)) { me=0;}
      
    // else if (count == nmol) { break;}
    // printf("count= %d\n ",counter);


              }
    

  }
printf("core %i count= %d\n ",me,counter);


////////////////// end fourth timer ////////////
  cputime= clock()-cputime;      // calculate  cpu clock time as difference of times after-before
  gettimeofday(&end, NULL);
  calc_wall_time = ((end.tv_sec  - start.tv_sec)+(end.tv_usec - start.tv_usec)/1e6);
  calc_cpu_time = (float) cputime/CLOCKS_PER_SEC;
////////////////// end fourth timer ////////////



// double Tserial10 = 0.891796;
// double Tserial100 = 61.748388;


Ttime = calc_wall_time  + read_wall_time + Bcast_wall_time;
cpu_total = Ttime*nproc;
theo_npair =  nmol*  (nmol-1)  /2;
// speedup = Tserial10 / Ttime;
// overhead = ( (Ttime*nproc/Tserial10) -1) *100 ;

MPI_Reduce(&calc_wall_time,&max_calc,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD); 
MPI_Reduce(&Ttime,&maxtime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD); 
MPI_Reduce(&Ttime,&mintime,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD); 
// MPI_Reduce(&cpu_total,&cpu_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); 

MPI_Allreduce(&energy,&total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 
MPI_Allreduce(&counter,&maxnpair,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
MPI_Allreduce(&counter,&minnpair,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
MPI_Allreduce(&counter,&actual_nparir,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); 


imbalance = 2.*(maxnpair-minnpair)/(0.+maxnpair+minnpair);


    


  if(me==0){
  printf("Total parallelized energy : %lf \n ",total);
  printf("Energy per molecule %lf \n",total/nmol);

  printf("before read wall time: %lf\n", before_read_wall_time);
  printf("before read CPU  time: %lf\n", before_read_cpu_time);
  printf("read wall time: %lf\n", read_wall_time);
  printf("read CPU  time: %lf\n", read_cpu_time);
  printf("read and Bcast wall time: %lf\n",  read_wall_time + Bcast_wall_time);
  printf("read and Bcast CPU  time: %lf\n",  read_cpu_time  + Bcast_cpu_time);
  printf("calc wall time: %lf\n", calc_wall_time);
  printf("calc CPU  time: %lf\n", calc_cpu_time);
  printf("Total wall  time: %lf\n", Ttime);
  // printf("Theoretical time: %lf\n", Tserial10/nproc);
  printf("MAX calc wall time: %lf\n", max_calc);

  printf("Theoretical npairs: %i\n", theo_npair);
  printf("Actual npairs: %i\n", actual_nparir);

  printf(" CPU Total time: %lf\n", cpu_total);



  printf("imbalance: %lf\n", imbalance);
  printf("max pairs: %i\n", maxnpair);
  printf("min pairs: %i\n", minnpair);
  printf("max read & calc  time: %lf\n", maxtime);
  printf("min read & calc time: %lf\n", mintime);
  // printf("speedup: %lf\n", speedup);
  // printf("overhead percentage: %lf\n", overhead);






 }

   fclose(fp);


MPI_Finalize(); // finalize MPI peacefully (the system would kill the processes otherwise)



}
