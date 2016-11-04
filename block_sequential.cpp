
#include<iostream>
#include<limits.h>
#include<stdio.h>
#include <ctime>
#include <omp.h>
using namespace std;

int main()
{
		int N,K,dim,hello;
		int cnt=0;
		clock_t begin,end;
		double elapsed_secs;
		cin>>N>>dim>>K;
		int i,k,j;
		float points[N][dim];
		float cluster[K][dim];
		float newCluster[K][dim];
		int membership[N],newClusterSize[K];
		float dist[K];
		int x,y;
		for(i=0;i<N;i++)
		{
			for(j=0;j<dim;j++)
				scanf("%f",&points[i][j]);
			membership[i]=-1;
		}
		for(k=0;k<K;k++)
		{
			for(j=0;j<dim;j++)
				cluster[k][j]=points[k][j];
			newClusterSize[k]=0;
			dist[k]=0;
		}
		int change,indx;
		float min_dist;
		clock_t time1=clock();
		double clustering_timing=omp_get_wtime();

		do
		{
			change=0;
			cnt++;
			for(i=0;i<N;i++)
			{
				indx=0;
				int min_dist=INT_MAX;
				for( k=0;k<K;k+=2)
				{
					for(j=0;j<dim;j+=2)
					{
					    for (x = k; x < min(k + 2, K); x++)
                            {
                                for (y = j; y < min(j + 2, dim); y++)
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
				newClusterSize[indx]++;
				for(int j=0;j<dim;j++)
				{
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
		printf("Total time taken for sequential= %lf\n",clustering_timing);
		for(k=0;k<K;k++)
		{
			for(j=0;j<dim;j++)
			{
                printf("%f ",cluster[k][j]);
			}
			cout<<endl;
		}
    }
