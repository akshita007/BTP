#include<iostream>
#include<limits.h>
#include<stdio.h>
#include <ctime>
#include <omp.h>
#include <unistd.h>   
#include <stdlib.h>
#include <sys/types.h>  /* open() */
#include <sys/stat.h>
#include <getopt.h>

using namespace std;
float points[10000000][10];

static void usage(char *argv0) {
    const char *help =
        "Usage: %s [switches] -i filename -n num_clusters\n"
        "       -n Number of points    : No of input data points\n"
        "       -d Dimensionality    : Dimensionality of each point\n"	
        "       -i filename    : file containing data to be clustered\n"
        "       -k num_clusters: number of clusters (K must > 1)\n"
        "       -h             : print this help information\n";
    fprintf(stderr, help, argv0);
    exit(-1);
}



int main(int argc, char *argv[])
{
		int N,K,dim;
		int option;
		char   *input_file;
		clock_t begin,end;
		double elapsed_secs;
		FILE *pFile;
		 while ((option = getopt(argc, argv,"n:d:k:i:")) != -1) {
        		switch (option) {
             			case 'n' : N=atoi(optarg);
				break;
             			case 'd' : dim=atoi(optarg);
                 		break;
             			case 'k' : K=atoi(optarg);
                 		break;
             			case 'i' : input_file=optarg;
                 		break;
             			default: usage(argv[0]); 
                 		exit(EXIT_FAILURE);
        		}
    		}
		pFile = fopen (input_file,"r");
		int i,k,j;
		if ( pFile == NULL) {
    			printf("\nCould not open file");
			return 0;
  		}
		float cluster[K][dim];
		float newCluster[K][dim];
		int membership[N],newClusterSize[K];
		float dist[K];
		int x,y;
		for(i=0;i<N;i++)
		{
			for(j=0;j<dim;j++)
				fscanf(pFile,"%f",&points[i][j]);
			membership[i]=-1;
		}
		fclose(pFile);
		for(k=0;k<K;k++)
		{
			for(j=0;j<dim;j++){
				cluster[k][j]=points[k][j];
				newCluster[k][j]=0;
			}
			newClusterSize[k]=0;
			dist[k]=0;
		}
		int change,indx;
		float min_dist;
		clock_t time1=clock();
		double clustering_timing=omp_get_wtime();
		int block_size=2;
		do
		{
			change=0;
			#pragma omp parallel for private(j,k,x,y,indx,dist) schedule(static) reduction(+:change)
			for(i=0;i<N;i++)
			{
				indx=0;
				int min_dist=INT_MAX;
				for( k=0;k<K;k+=block_size)
				{
					for(j=0;j<dim;j+=block_size)
					{
					    for (x = k; x < min(k + block_size, K); x++)
                            {
                                for (y = j; y < min(j + block_size, dim); y++)
                                {
                                    dist[x]+=(cluster[x][y]-points[i][y])*(cluster[x][y]-points[i][y]);
                                    if(j==dim-1)
                                    {
                                        if(dist[x]<min_dist)
                                        {
                                            min_dist=dist[x];
                                            indx=x;
                                        }
                                        dist[x]=0;
                                    }
                                }
                            }
					}
				}//calculated distance with every point
				if(membership[i]!=indx)
                    change++;

				membership[i]=indx;
				#pragma omp atomic
				newClusterSize[indx]++;
				for(int j=0;j<dim;j++)
				{
				   #pragma omp atomic
                    newCluster[indx][j]+=points[i][j];
				}
            }
			for(k=0;k<K;k++)
			{
				for(int j=0;j<dim;j++)
				{
                    cluster[k][j]=newCluster[k][j]/newClusterSize[k];
                    newCluster[k][j]=0;
				}
				newClusterSize[k]=0;
			}
			//end=clock();
			//elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			//printf("Recalculation stage time = %lf\n",elapsed_secs);
		}while(change>0);
		clustering_timing=omp_get_wtime()-clustering_timing;
		printf("Total time taken for blocking= %lf\n",clustering_timing);
		for(k=0;k<K;k++)
		{
			for(j=0;j<dim;j++)
			{
                printf("%f ",cluster[k][j]);
			}
			cout<<endl;
		}
    }
