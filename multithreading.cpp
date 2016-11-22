
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
		pFile = fopen ("output.txt","r");
		cin>>N>>dim>>K;
		int i,k,j;
		float cluster[K][dim];
		float newCluster[K][dim];
		int membership[N],newClusterSize[K];
		//#pragma omp for
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
		}
		int change,indx;
		float min_dist,dist;
		clock_t time1=clock();
		double cluster_timing;
		cluster_timing = omp_get_wtime();
		do
		{
			change=0;
			cnt++;
			begin = clock();
			#pragma omp parallel for private(j,k,indx,dist) schedule(static) reduction(+:change) shared(points,cluster,membership,newCluster,newClusterSize)
			for(i=0;i<N;i++)
			{
				indx=0;
				double  min_dist=INT_MAX;
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
				#pragma omp atomic
				newClusterSize[indx]++;
				for(int j=0;j<dim;j++){
					#pragma omp atomic
					newCluster[indx][j]+=points[i][j];
				}
			}
			for(k=0;k<K;k++)
			{
				for(int j=0;j<dim;j++){
						cluster[k][j]=newCluster[k][j]/newClusterSize[k];
				newCluster[k][j]=0;
				}
				newClusterSize[k]=0;
			}
			//end=clock();
			//elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			//printf("Recalculation stage time = %lf\n",elapsed_secs);
		}while(change>0);
		clock_t time2=clock();
		elapsed_secs = double(time2-time1) / CLOCKS_PER_SEC;
		cluster_timing=omp_get_wtime()-cluster_timing;
		//printf("Total time taken for threaded= %lf\n",elapsed_secs);
		printf("Total time taken for threaded= %lf\n",cluster_timing);

		for(k=0;k<K;k++)
		{
			for(j=0;j<dim;j++){
                printf("%f ",cluster[k][j]);
			}
			cout<<endl;
		}


}

//calculated the
