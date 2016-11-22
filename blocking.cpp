#include<iostream>
#include<limits.h>
#include<stdio.h>
#include <ctime>
#include <omp.h>
using namespace std;
float points[10000000][10];
int main()
{
		int N,K,dim,hello;
		int cnt=0;
		clock_t begin,end;
		double elapsed_secs;
		FILE *pFile;
		pFile = fopen ("Large_Input","r");
		cin>>N>>dim>>K;
		int i,k,j;
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
			cnt++;
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
