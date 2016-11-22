

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
#include<float.h>

using namespace std;
float points[10][3000000];
int membership[3000000];
int newmembership[3000000];
float min_dist[3000000],dist[3000000];

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
		if ( pFile == NULL) {
    			printf("\nCould not open file");
			return 0;
  		}

		int i,k,j;
		float cluster[dim][K];
		float newCluster[dim][K];
		int newClusterSize[K];
		for(i=0;i<N;i++)
		{
			for(j=0;j<dim;j++)
				fscanf(pFile,"%f",&points[j][i]);
			membership[i]=-1;
			//cout<<"Input "<<i<<endl;
		}
		fclose(pFile);
		for(j=0;j<dim;j++)
		{
			for(k=0;k<K;k++)
				cluster[j][k]=points[j][k];
		}
		int change=0,indx;
		double clustering_timing=omp_get_wtime();
		do
		{
			change=0;
			for(i=0;i<N;i++)
			{
				min_dist[i]=FLT_MAX;
				dist[i]=0;
			}
			for(k=0;k<K;k++)
			{
				//for each cluster compute it distance with each point
				newClusterSize[k]=0;
				for(j=0;j<dim;j++)
				{
					newCluster[j][k]=0;
					for(i=0;i<N;i++)
					{
						dist[i]+=(cluster[j][k]-points[j][i])*(cluster[j][k]-points[j][i]);
					}
				}

				for(i=0;i<N;i++)
				{
					if(min_dist[i]>dist[i])//The point i belongs to the cluster k
					{
						newmembership[i]=k;
						min_dist[i]=dist[i];
					}
					dist[i]=0;
				}
			}
			for(j=0;j<dim;j++){
					for(i=0;i<N;i++){
							if(j==0){
						newClusterSize[newmembership[i]]+=1;
						if(newmembership[i]!=membership[i]){
							change++;
						membership[i]=newmembership[i];}
							}
						newCluster[j][membership[i]]+=points[j][i];

					}
			}
			for(k=0;k<K;k++)
			{
				for(j=0;j<dim;j++){
					cluster[j][k]=newCluster[j][k]/newClusterSize[k];
				}
			}
		}while(change>0);
		clustering_timing=omp_get_wtime()-clustering_timing;
		printf("Total time taken for centroid centric= %lf\n",clustering_timing);
		for(k=0;k<K;k++)
		{
			for(j=0;j<dim;j++){
                printf("%f ",cluster[j][k]);
			}
			cout<<endl;

		}
}
