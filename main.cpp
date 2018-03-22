#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "stdio.h"
#include <ios>

#define MASTER 0              
#define FROM_MASTER 1     
#define FROM_WORKER 2      

using namespace std;

float checkforrepeats(int **sdna,int **x,int dim,int w,int go_on){
    for(int i=0;i<w;i++){
        int cnt=0;
        for(int j=0;j<dim;j++)
          if(x[0][j]==sdna[i][j])
             cnt++;
          if(cnt == dim){
            go_on=1;
            break;
        }
    }
    return go_on;
}

void dna_trajectory(int min,int max, int length_dna,int **sdna,int dim) {
    int ** x;
    x = new int *[length_dna+1];
    for (int i=0; i<length_dna+1; i++) {
		x[i] = new int [dim];
    }
	
    for (int i=0; i<length_dna+1; i++) {
        for(int j=0;j<dim;j++){
            x[i][j] = 0;
        }
    }
	
    int w=1,no_move=0;
    while( w < length_dna) {
		for(int i=0;i<dim;i++){
			float rval = (rand()%10)+1;
			
			if(rval < 5){
				if(static_cast<int>(sdna[w-1][i]-1) >= min)
					x[0][i] = static_cast<int>(sdna[w-1][i] - 1);
				
				else
					x[0][i] = static_cast<int>(sdna[w-1][i] + 1);
			}
			else {
				if(static_cast<int>(sdna[w-1][i]+1) <= max)
					x[0][i] =static_cast<int>(sdna[w-1][i]+1);
				else
					x[0][i] = static_cast<int>(sdna[w-1][i]-1);
			}
			
		}
		int go_on=0;
		
		go_on=checkforrepeats(sdna,x,dim,w,go_on);
		
		if(go_on==0){
			for(int k=0;k<dim;k++)
				sdna[w][k]=x[0][k];
			w++;
		}
		else {
			no_move++;
		}
		
		if(no_move > 10) {
			no_move=0;
			w=1;
		}
    }
}

int main (int argc, char *argv[])
{
	int num_avg=100000,dim=2,min=0,max=500;
	int bp_start=5,bp_end=200;
	int bps=bp_end-bp_start+1;
	int	numtasks,              /* number of tasks in partition */
	taskid,                /* a task identifier */
	numworkers,            /* number of worker tasks */
	source,                /* task id of message source */
	dest,                  /* task id of message destination */
	mtype,                 /* message type */
	rows,                  /* rows of matrix A sent to each worker */
	averow, extra, offset, /* used to determine rows sent to each worker */
	i, j, k, rc;           /* misc */
	int nm;
	
	
	MPI_Status status;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	
	numworkers = numtasks-1;
	
	char file[20];
	sprintf(file,"SAW_%d.txt",dim);
	ofstream out(file);
	
	srand(time(NULL));
	
	for(int l=0;l<bps;l++){
		double	a[num_avg][1], c[num_avg][2];
		int length_dna=l+bp_start;
		int ** sdna_one , ** sdna_two;
		sdna_one = new int *[length_dna];
		sdna_two = new int *[length_dna];
		for (int k=0; k<length_dna; k++) {
			sdna_one[k] = new int [dim];
			sdna_two[k] = new int [dim];
		}
		
	 /**************************** master task ************************************/
		if (taskid == MASTER)
		{
			for (i=0; i<num_avg; i++)
				a[i][0]= i;
			
			/* Send matrix data to the worker tasks */
			averow = num_avg/numworkers;
			extra = num_avg%numworkers;
			offset = 0;
			mtype = FROM_MASTER;
			for (dest=1; dest<=numworkers; dest++)
			{
				rows = (dest <= extra) ? averow+1 : averow;   
				MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
				MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
				MPI_Send(&a[offset][0], rows, MPI_DOUBLE, dest, mtype,MPI_COMM_WORLD);
				
				offset = offset + rows;
			}
			
			/* Receive results from worker tasks */
			mtype = FROM_WORKER;
			for (i=1; i<=numworkers; i++)
			{
				source = i;
				MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
				MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
				MPI_Recv(&c[offset][0], rows*2, MPI_DOUBLE, source, mtype,MPI_COMM_WORLD, &status);

			}
			
			/* Print results */
			int cont=0,corr_cont=0;
			for (i=0; i<num_avg; i++)
			{
				cont=cont+c[i][0];
				corr_cont=corr_cont+c[i][1];
			}
			float no_cont=cont/(num_avg+0.0);
			float no_corr_cont=corr_cont/(num_avg+0.0);
			float prob = static_cast<float>(no_corr_cont/no_cont);
			
			out << length_dna << "," << no_cont << "," << no_corr_cont << "," << prob << std::endl;
		}
		
		
		/**************************** worker task ************************************/
		if (taskid > MASTER)
		{
			mtype = FROM_MASTER;
			MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&a, rows, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
			
			for(int n=0;n<rows;n++){
				for (int k=0; k<length_dna; k++) {
					for(int p=0;p<dim;p++){
						sdna_one[k][p] = 0;
						sdna_two[k][p] = 0;
					}
				}

				int tot_cont=0, tot_corr_cont=0;
				dna_trajectory(min,max,length_dna,sdna_one,dim);
				dna_trajectory(min,max,length_dna,sdna_two,dim);
				
				for(i=0;i<length_dna;i++){
				//	std::cout << sdna_one[i][0] << ","<<sdna_one[i][1] << std::endl;
				//	std::cout << sdna_two[i][0] << ","<<sdna_two[i][1] << std::endl;
					for(j=0;j<length_dna;j++){
						int cnt=0;
						for(int k=0;k<dim;k++)
							if(sdna_one[i][k]==sdna_two[j][k])
								cnt++;
						if(cnt==dim)
							tot_cont++;
						if(cnt==dim && i==j)
							tot_corr_cont++;
					}
				}
				c[n][0]=tot_cont;
				c[n][1]=tot_corr_cont;
				
			}
			mtype = FROM_WORKER;
			MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
			MPI_Send(&c, rows*2, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
		}
	}
	MPI_Finalize();
}
