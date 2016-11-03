

#include<iostream>
#include<limits.h>
#include<float.h>
#include<omp.h>
#include<stdio.h>
using namespace std;

int main()
{
        int N;
		int K,dim;
		//cout<<"Enter the number of points, dimensionality of each and the number of clusters";
		//cin>>N>>
		cin>>N>>dim>>K;
		int i,k,j;
		float points[dim][N];//the objects will be stored dimension wise
		float cluster[dim][K];
		float newCluster[dim][K];
		omp_lock_t writelock[N];
        for(i=0;i<N;i++)
            omp_init_lock(&writelock[i]);

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
			change=0;
			for(i=0;i<N;i++)
			{
				min_dist[i]=FLT_MAX;
				dist[i]=0;
			}
			#pragma omp parallel for private(dist,i,j)
			for(k=0;k<K;k++)
			{
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
				//#pragma omp critical
				//{
				omp_set_lock(&(writelock[i]));
					if(min_dist[i]>dist[i])
					{
						newmembership[i]=k;
						min_dist[i]=dist[i];
					}
                omp_unset_lock(&(writelock[i]));

                //}
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
