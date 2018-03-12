#include <stdio.h>
#include <string.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>



//---------------- initinalisation and declaration-------------
int Size = 2000000;
int comm_sz, rank; 
double t5, t6;


void SerialAlgorithm();
void ParallelAlgorithm();


//--------------------main-----------------------

int main(void) {
    
    //-----init mpi---
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
     if (rank == 0) {
        SerialAlgorithm();
    }
     if (comm_sz > 1) {
        
       ParallelAlgorithm();
      
      }
        else {
        printf("\n\n*** Please choose 2 or more processes.\n\n");
    }

MPI_Finalize();
return 0;
    
}

//-----------Serial Algorithm starts here----------------------------
void SerialAlgorithm(){

printf("\n\n********************Start the serial Algorithm********************\n\n");

//---------------start the clock----------
double t1, t2;
//MPI_Init(0,0);

t1 =  MPI_Wtime();

//---------------Declarations---------------------- 

double* a = malloc((Size) * sizeof(double));
double* b = malloc((Size) * sizeof(double));
double suma = 0;    
double sumb = 0;
double stdeva = 0;  
double stdevb = 0;
double arraya = 0; 
double arrayb = 0;
double arraysab = 0;
double pearson = 0; 

//---------------Populate the arrays-------------
 
  for(int i = 0; i < Size; i++){

  
    a[i]= sin(i);
    suma += a[i];
    
    b[i]= sin(i+5);
    sumb += b[i];
  }
//---------------Calculate the mean ----------------

double meana = suma / Size;
double meanb = sumb / Size;

printf("\n**********Mean of arrays a and b in serial code**********\n\n");

printf("*** The mean of array a is: %.10f\n", meana);
printf("*** The mean of array b is: %.10f\n\n", meanb);

/*---------------Calculate the numerator for both stadard deviation and pearson 
                    in order to avoid repeating the same function twice-----------*/


  for (int j = 0; j < Size; j++){
        
        arraya += (a[j] - meana) * (a[j] - meana);
        arrayb += (b[j] - meanb) * (b[j] - meanb);
        arraysab += (a[j] - meana) * (b[j] - meanb);   
  }

//---------------standard devaition ---------------

printf("**********Standard deviation of arrays a and b in serial code:**********");

stdeva = sqrt(arraya / Size);
printf("\n\n*** The standard deviation of array a is: %.10f\n", stdeva);

stdevb = sqrt(arrayb / Size);
printf("*** The standard deviation of array b is: %.10f\n\n\n", stdevb);

//--------------------Pearson-----------------------

printf("**********Pearson's correlation of arrays a and b in serial code**********");

arraysab /= Size;
pearson = arraysab / (stdeva * stdevb);
printf("\n\n*** Pearson's correlation coefficient is: %.10f", pearson);

//---------------Measure the time----------- 
t2 = MPI_Wtime();
t5 = t2-t1;
printf("\n\nTime measured: %.2fs\n\n", t5);fflush(stdout);

printf("********************End of Serial Algorithm********************\n");

free(a);
free(b);
}



//---------------Parallel algorithm--------
void ParallelAlgorithm(){
double t3, t4;

  if (rank == 0){
    
     t3 =  MPI_Wtime();
     printf("\n\n\n\n********************Start the parallel algorithm********************\n\n");
  
  }

//---------------Declarations---------------------- 

double suma = 0;
double sumb = 0;
double meana = 0;
double meanb = 0;
double stdeva = 0;
double stdevb = 0;
double arraya = 0;
double arrayb = 0;
double arraysab = 0;
double pearson = 0; 



//------------Check if the nb of element can be devided by the number of processes---
// If yes, run the program normally 
//If no, some process will do more work than others, this is done by adding more elements to some process
//More specifically, the processes that have a rank inferior to the rest will work more then others.
 

int Sizebis = 0; // New Size depending on the current process 
int rest = Size % comm_sz;
    
    if(rest == 0){
      Sizebis = Size/comm_sz;  // no additional work for these processes 
      }
  else{
      Sizebis = Size/comm_sz; 
       if (rank < rest){
            Sizebis++;    // additional work for these processes.
        }
    }

//---------------declaring the arrays where the elements will be stored---
//---------------the size of the arrays depend on whther it's a process that needs additional work or not---

double* abis = malloc(Sizebis * sizeof(double));
double* bbis = malloc(Sizebis * sizeof(double));

//---- to avoid making the same calculation for each process, i calculate the offset for each process 
//---- in this case, instead of restarting from zero a process will start from a specfic element in the table 

 double offset = 0;
    if (rank < rest) {
        offset = rank * Sizebis;
    }
    else {
        offset = (rank * Sizebis) + rest;
    }

//---------------Populate the arrays-------------
 
  for(int i = 0; i < Sizebis; i++){
  
    abis[i]= sin(i + offset);
    suma += abis[i];
    
    bbis[i]= sin(i+5+offset);
    sumb += bbis[i];
  }
  
// send it to the root 
double sumabis = 0;
double sumbbis = 0;

  MPI_Reduce(&suma, &sumabis, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sumb, &sumbbis, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  
//--------------Calculate the mean ----------------

if (rank == 0) {
             
        meana = sumabis / Size;
        meanb = sumbbis / Size;

        printf("\n**********Mean of arrays a and b in Parallel code**********\n\n");

        printf("*** The mean of array a is: %.10f\n", meana);
        printf("*** The mean of array b is: %.10f\n\n", meanb);
    }


//------------- Send the mean calculated to all the proc-------

    MPI_Bcast(&meana, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&meanb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    

/*---------------Calculate the numerator for both stadard deviation and pearson 
                    in order to avoid repeating the same function twice-----------*/

  for (int j = 0; j < Sizebis; j++){
        
        arraya += (abis[j] - meana) * (abis[j] - meana);
        arrayb += (bbis[j] - meanb) * (bbis[j] - meanb);
        arraysab += (abis[j] - meana) * (bbis[j] - meanb);   
  }

double numa = 0;
double numb = 0;
double numab = 0;

    MPI_Reduce(&arraya, &numa, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&arrayb, &numb, 1, MPI_DOUBLE, MPI_SUM, 1, MPI_COMM_WORLD);
    MPI_Reduce(&arraysab, &numab, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


free(bbis);
free(abis);


//---------------standard devaition ---------------

 
    if (rank == 1) {
        
        stdevb = sqrt(numb / Size);
        MPI_Send(&stdevb, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    if (rank == 0) {
    
        printf("*********Standard deviation of arrays a and b in serial code**********");
      
        double stdeva = sqrt(numa / Size);
        printf("\n\n*** The standard deviation of array a is: %.10f\n", stdeva);
        
        MPI_Recv(&stdevb, 1, MPI_DOUBLE, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("*** The standard deviation of array b is: %.10f\n\n\n", stdevb);

//------------------Pearson-----------------------

    printf("**********Pearson's correlation of arrays a and b in serial code**********");
         
          numab /= Size;
          pearson = numab / (stdeva * stdevb);
          printf("\n\n*Pearson's correlation coefficient is: %.10f", pearson);


//--------------Measure the time--------------------- 
t4 = MPI_Wtime();
t6 = t4 - t3;
printf("\n\nTime measured: %1.2fs\n\n", t6);fflush(stdout);
printf("***********End of the parallel Code************\n");
  
  printf("\n*** Seedup: %.10fs \n\n", (t5 - t6));
  printf("\n\n*** Percentage Speedup: %.0f%%\n", (1-(t6/t5))*100);

  }
}
