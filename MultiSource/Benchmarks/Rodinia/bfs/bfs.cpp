/*
- Modified by Pankaj Kukreja
- Indian Institute of Technology Hyderabad, India
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
//Structure to hold a node information

// Gap in which array should be printed
#define GAP 20

struct Node
{
	int source;
	int edge_id;
};
void BFSGraph(char *filename);

////////////////////////////////////////////////////////////////////////////////
// Main Program
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[]) 
{

	if(argc!=2)
	{
		fprintf(stderr,"Usage: %s <input_file>\n", argv[0]);
		exit(0);
	}
	BFSGraph(argv[1]);
	return 0;
}

void BFSGraph(char* filename)
{
	FILE *fp;
    int no_of_nodes = 0;
    int edge_list_size = 0;

	//Read in Graph from a file
	fp = fopen(filename,"r");
	if(!fp)
	{
		printf("Error Reading graph file\n");
		exit(0);
	}
	int source = 0;
	fscanf(fp,"%d",&no_of_nodes);
    
	// allocate host memory
	Node * __restrict__ h_graph_nodes = (Node*) malloc(sizeof(Node)*no_of_nodes);
	bool * __restrict__ h_graph_mask = (bool*) malloc(sizeof(bool)*no_of_nodes);
	bool * __restrict__ h_updating_graph_mask = (bool*) malloc(sizeof(bool)*no_of_nodes);
	bool * __restrict__ h_graph_visited = (bool*) malloc(sizeof(bool)*no_of_nodes);
	int start, edgeid;   
	// initalize the memory
	for(int i = 0; i < no_of_nodes; i++) 
	{
		fscanf(fp,"%d %d",&start,&edgeid);
		h_graph_nodes[i].source = start;
		h_graph_nodes[i].edge_id = edgeid;
		h_graph_mask[i]=false;
		h_updating_graph_mask[i]=false;
		h_graph_visited[i]=false;
	}
	//read the source node from the file
	fscanf(fp,"%d",&source);
	// source=0; //tesing code line
	//set the source node as true in the mask
	h_graph_mask[source]=true;
	h_graph_visited[source]=true;
	fscanf(fp,"%d",&edge_list_size);
	int id,cost;
	int* __restrict__ h_graph_edges = (int*) malloc(sizeof(int)*edge_list_size);
	
	//  Here we changed a bit
	for(int i=0; i < edge_list_size ; i++)
	{
		fscanf(fp,"%d",&id);
		fscanf(fp,"%d",&cost);
		h_graph_edges[i] = id;
	}

	if(fp)
		fclose(fp);    

	// allocate mem for the result on host side
	int* __restrict__ h_cost = (int*) malloc(sizeof(int)*no_of_nodes);
	for(int i=0;i<no_of_nodes;i++)
		h_cost[i]=-1;
	h_cost[source]=0;

	int k=0;
	bool stop;
	do
    {
        stop=false;
        for(int _k = 0; _k < no_of_nodes; _k++ )
        {
            if(h_graph_mask[_k]==true)
            { 
                h_graph_mask[_k]=false;
                for(int i=h_graph_nodes[_k].source;
                	i<(h_graph_nodes[_k].edge_id + h_graph_nodes[_k].source); 
                	i++)
                {
                    int id = h_graph_edges[i];
                    if(!h_graph_visited[id])
                    {
                        h_cost[id]=h_cost[_k]+1;
                        h_updating_graph_mask[id]=true;
                    }
                }
            }
        }
        for(int _k=0; _k< no_of_nodes ; _k++ )
        {
            if (h_updating_graph_mask[_k] == true)
            {
                h_graph_mask[_k]=true;
                h_graph_visited[_k]=true;
                stop=true;
                h_updating_graph_mask[_k]=false;
            }
        }
        k++;
    }
	while(stop);
	//Store the result into a file
	for(int i=0;i<no_of_nodes;i++)
	{
		if(i % GAP == 0)
		{
			fprintf(stdout,"%d %d\n",i,h_cost[i]);
		}
	}
	
	// cleanup memory
	free( h_graph_nodes);
	free( h_graph_edges);
	free( h_graph_mask);
	free( h_updating_graph_mask);
	free( h_graph_visited);
	free( h_cost);
}
