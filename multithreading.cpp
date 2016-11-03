
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
		//cout<<"Enter the number of points, dimensionality of each and the number of clusters";
		cin>>N>>dim>>K;
		int i,k,j;
		float points[N][dim];
		float cluster[K][dim];
		float newCluster[K][dim];
		int membership[N],newClusterSize[K];
		//#pragma omp for
		for(i=0;i<N;i++)
		{
			//scanf("%d",&j);//ignoring the first number
			for(j=0;j<dim;j++)
				scanf("%f",&points[i][j]);
			membership[i]=-1;
		}
		//#pragma omp for
		for(k=0;k<K;k++)
		{
			for(j=0;j<dim;j++)
				cluster[k][j]=points[k][j];
			newClusterSize[k]=0;
		}
		int change,indx;
		float min_dist,dist;
		clock_t time1=clock();
		double cluster_timing;
		cluster_timing = omp_get_wtime();
		do
		{
			//objects assignment

			change=0;
			cnt++;
			begin = clock();
			#pragma omp parallel for private(i,j,k,indx,dist) schedule(static) reduction(+:change)
			for(i=0;i<N;i++)
			{
				indx=0;
				//int tid = omp_get_thread_num();
				int min_dist=INT_MAX;
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
				//printf("Thread id=%d\n",tid);
			}
			//end = clock();
			//elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			//printf("Assignment stage time = %lf count=%d\n",elapsed_secs,cnt);
			//begin=clock();
			//#pragma omp parallel for
			for(k=0;k<K;k++)
			{
				for(int j=0;j<dim;j++){
						cluster[k][j]=newCluster[k][j]/newClusterSize[k];
				newCluster[k][j]=0;
				}
				//cout<<newClusterSize[k]<<endl;
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
