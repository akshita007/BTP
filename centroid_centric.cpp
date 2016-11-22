

#include<iostream>
#include<limits.h>
#include<float.h>
#include<omp.h>
#include<stdio.h>
using namespace std;
float points[10][3000000];
int membership[3000000];
int newmembership[3000000];
float min_dist[3000000],dist[3000000];
int main()
{
		int N,K,dim;
		//cout<<"Enter the number of points, dimensionality of each and the number of clusters";
		cin>>N>>dim>>K;
		int i,k,j;
		float cluster[dim][K];
		float newCluster[dim][K];
		int newClusterSize[K];
		FILE *pFile;
		pFile = fopen ("output.txt","r");
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
                cout<<cluster[j][k]<<" ";
			}
			cout<<endl;

		}
}
