#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <time.h>



//---------------- initinalisation and declaration-------------
int Size = 5000000;
double t5, t6;


void SerialAlgorithm();
void ParallelAlgorithm();


//--------------------main-----------------------

int main(void) {

        SerialAlgorithm();
        ParallelAlgorithm();

return 0; 
}

//-----------Serial Algorithm starts here----------------------------
void SerialAlgorithm(){

printf("\n\n********************Start the serial Algorithm********************\n\n");

//---------------start the clock----------
double t1, t2;

t1 =  omp_get_wtime();

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
t2 = omp_get_wtime();
t5 = t2-t1;
printf("\n\nTime measured: %.2fs\n\n", t5);fflush(stdout);

printf("********************End of Serial Algorithm********************\n");

free(a);
free(b);
}



//----------------------Parallel algorithm------------------
void ParallelAlgorithm(){

printf("\n\n********************Start the Parallel Algorithm********************\n\n");

//---------------start the clock----------
double t3, t4;

t3 =  omp_get_wtime();
int thds = 1;
//---------------Declarations---------------------- 

double* a = malloc((Size) * sizeof(double));
double* b = malloc((Size) * sizeof(double));
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

//---------------Populate the arrays-------------
 
#pragma omp parallel for schedule(static) num_threads(thds) reduction(+:suma) reduction(+:sumb) 

  for(int i = 0; i < Size; i++){
  
    a[i]= sin(i);
    suma += a[i];
    
    b[i]= sin(i+5);
    sumb += b[i];
  }
//---------------Calculate the mean ----------------

printf("\n**********Mean of arrays a and b in parallel code**********\n\n");

#pragma omp parallel sections num_threads(thds)
{
#pragma omp section
{
double meana = suma / Size;
printf("*** The mean of array a is: %.10f\n", meana);
}
#pragma omp section
{
double meanb = sumb / Size;
printf("*** The mean of array b is: %.10f\n\n", meanb);
  }
}

/*---------------Calculate the numerator for both stadard deviation and pearson 
                    in order to avoid repeating the same function twice-----------*/

#pragma omp parallel for schedule(static) num_threads(thds) reduction(+:arraya) reduction(+:arrayb) reduction(+:arraysab)

  for (int j = 0; j < Size; j++){
        
        arraya += (a[j] - meana) * (a[j] - meana);
        arrayb += (b[j] - meanb) * (b[j] - meanb);
        arraysab += (a[j] - meana) * (b[j] - meanb);   
  }

//---------------standard devaition ---------------

printf("**Standard deviation of arrays a and b in parallel code:**\n\n");

#pragma omp parallel sections num_threads(thds)
{
#pragma omp section
{
stdeva = sqrt(arraya / Size);
printf("*** The standard deviation of array a is: %.10f\n\n", stdeva);
}
#pragma omp section
{
stdevb = sqrt(arrayb / Size);
printf("*** The standard deviation of array b is: %.10f\n\n", stdevb);
  }
}
//--------------------Pearson-----------------------

printf("**********Pearson's correlation of arrays a and b in parallel code**********");

arraysab /= Size;
pearson = arraysab / (stdeva * stdevb);
printf("\n\n*** Pearson's correlation coefficient is: %.10f", pearson);

//---------------Measure the time----------- 
t4 = omp_get_wtime();
t6 = t4 - t3;
printf("\n\n Time measured: %1.2fs\n\n", t6);fflush(stdout);

printf("********************End of Parallel Algorithm********************\n");  
  printf("\n*** Seedup: %.10fs \n", (t5 - t6));
  printf("*** Percentage Speedup: %.0f%%\n\n", (1-(t6/t5))*100);



free(a);
free(b);
}
