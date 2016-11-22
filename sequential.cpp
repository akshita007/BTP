#include<iostream>
#include<limits.h>
#include<stdio.h>
#include <ctime>
#include <omp.h>
using namespace std;
double points[10000000][10];
int main()
{
		int N,K,dim,hello;
		int cnt=0;
		clock_t begin,end;
		double elapsed_secs;
		FILE *pFile;
		pFile = fopen ("output.txt","r");
		cin>>N>>dim>>K;
		int i,k,j;
		//float points[N][dim];
		double cluster[K][dim];
		double newCluster[K][dim];
		int membership[N],newClusterSize[K];
		for(i=0;i<N;i++)
		{
			//scanf("%d",&j);//ignoring the first number
			for(j=0;j<dim;j++)
				fscanf(pFile,"%lf",&points[i][j]);
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
		}
		int change,indx;
		double min_dist,dist;
		clock_t time1=clock();
		double clustering_timing=omp_get_wtime();
		do
		{
			change=0;
			cnt++;
			for(i=0;i<N;i++)
			{
				indx=0;
				min_dist=INT_MAX;
				for( k=0;k<K;k++)
				{
					dist=0;
					for(j=0;j<dim;j++)
					{
						dist+=(cluster[k][j]-points[i][j])*(cluster[k][j]-points[i][j]);
					}
					if(dist<min_dist)
					{
						min_dist=dist;
						indx=k;
					}
				}
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
