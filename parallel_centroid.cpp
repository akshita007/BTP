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
float min_dist[3000000];
//float dist[3000000];

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
		float cluster[dim][K];
		float newCluster[dim][K];
float dist[N];
		omp_lock_t writelock[N];
        for(i=0;i<N;i++)
            omp_init_lock(&writelock[i]);
 
		int newClusterSize[K];
		for(i=0;i<N;i++)
		{
			//cin>>j;//ignoring the first number
			for(j=0;j<dim;j++)
				fscanf(pFile,"%f",&points[j][i]);
			membership[i]=-1;
		}//random initialization
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
			#pragma omp parallel for private(dist,i,j,k) shared(min_dist,writelock,newClusterSize,newCluster) schedule(static)
			for(k=0;k<K;k++)
			{
				//cout<<"Implemented"<<k<<endl;
				newClusterSize[k]=0;
				//cout<<"Hey";
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
				omp_set_lock(&(writelock[i]));
					if(min_dist[i]>dist[i])
					{
						newmembership[i]=k;
						min_dist[i]=dist[i];
					}
                		omp_unset_lock(&(writelock[i]));
				dist[i]=0;
				}
			}
			for(j=0;j<dim;j++)
			{
                for(i=0;i<N;i++)
                {
                    if(j==0)
                    {
                        newClusterSize[newmembership[i]]+=1;
						if(newmembership[i]!=membership[i])
						{
                            change++;
                            membership[i]=newmembership[i];
                        }
                    }
                    newCluster[j][membership[i]]+=points[j][i];
                }
			}
			#pragma omp parallel for private(j)
			for(k=0;k<K;k++)
			{
				for(j=0;j<dim;j++)
					cluster[j][k]=newCluster[j][k]/newClusterSize[k];
			}
	//cout<<"Implemented"<<endl;
        }while(change>0);
		clustering_timing=omp_get_wtime()-clustering_timing;
		printf("Total time taken for centroid centric= %lf\n",clustering_timing);
		for(k=0;k<K;k++)
		{
            for(j=0;j<dim;j++)
                cout<<cluster[j][k]<<" ";
			cout<<endl;
        }
}
