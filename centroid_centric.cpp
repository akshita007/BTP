

#include<iostream>
#include<limits.h>
#include<float.h>
#include<omp.h>
#include<stdio.h>
using namespace std;

int main()
{
		int N,K,dim;
		//cout<<"Enter the number of points, dimensionality of each and the number of clusters";
		cin>>N>>dim>>K;
		int i,k,j;
		float points[dim][N];//the objects will be stored dimension wise
		float cluster[dim][K];
		float newCluster[dim][K];
		int membership[N],newClusterSize[K],newmembership[N];
		for(i=0;i<N;i++)
		{
			//cin>>j;//ignoring the first number
			for(j=0;j<dim;j++)
				cin>>points[j][i];
			membership[i]=-1;
		}//random initialization
		for(j=0;j<dim;j++)
		{
			for(k=0;k<K;k++)
				cluster[j][k]=points[j][k];
		}
		int change=0,indx;
		float min_dist[N],dist[N];
		double clustering_timing=omp_get_wtime();
		do
		{
			//objects assignment
			//set min dist
			//cout<<change<<endl;
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
					/*if(k==0)
						cout<<cluster[j][1]<<" ";*/
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
				//cout<<newClusterSize[k]<<" "<<endl;;
				for(j=0;j<dim;j++)
					cluster[j][k]=newCluster[j][k]/newClusterSize[k];
			}
			//now taking average
			//cout<<endl<<change<<endl;
		}while(change>0);
		clustering_timing=omp_get_wtime()-clustering_timing;
		printf("Total time taken for centroid centric= %lf\n",clustering_timing);
		for(k=0;k<K;k++)
		{
			//cout<<newClusterSize[k]<<" "<<endl;;
			for(j=0;j<dim;j++){
                cout<<cluster[j][k]<<" ";
			}
			cout<<endl;

		}
}
