#include "gm_graph.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <map>
#include <vector>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <fstream>
#include "mmio.h"

using namespace std;

typedef int node_t;
typedef int edge_t;

struct timeval T1, T2;

vector< pair<node_t, node_t> > parse_adjacency_file(char* filename, bool undirected){
  ifstream adjFile(filename);
  vector< pair<node_t, node_t> > adjVector;
  pair<node_t, node_t> temp;
  int source, sink;
  while( adjFile >> source >> sink)
  {
    temp.first = (node_t) source;
    temp.second = (node_t) sink;
    adjVector.push_back(temp);
    if(undirected){
      temp.second = (node_t) source;
      temp.first = (node_t) sink;
      adjVector.push_back(temp);
    }
  }
  return adjVector;
}

node_t max_node(std::vector<pair<node_t, node_t> > list){
  node_t mNode = 0;
  for(std::vector< pair<node_t, node_t> >::iterator it = list.begin(); it!= list.end(); ++it)
  {
    if((*it).first > mNode)
    {
      mNode = (*it).first;
    }
    if((*it).second > mNode)
    {
      mNode = (*it).first;
    }
  }
  return mNode;
}

void write_green_marl_file(char* filename, vector< pair<node_t, node_t> > edges, node_t N, edge_t M){

  // allocate space for degree counts
  // node_t* src = new node_t[M];
  // node_t* dst = new node_t[M];
  edge_t* deg = new edge_t[N];
  memset(deg, 0, sizeof(edge_t) * N);
  
  gm_graph* g = new gm_graph();
  g->prepare_external_creation(N, M);
  gettimeofday(&T1, NULL);

  printf("Calculating degrees.\n");
  for(unsigned int i = 0; i<M; ++i) //iterate over edges - should use iterator :/
  {
    edges[i].first = edges[i].first;
    edges[i].second = edges[i].second;
    deg[edges[i].first]++;
  }

  printf("Setting gm_graph degrees\n");
  //manually manipulate the sparse internal graph format
  //see graph_gen.cc in $GREEN_MARL/apps/output_cpp/gm_graph/src
  g->begin[0] = 0;
  for (node_t i = 1; i <= N; i++) {
    g->begin[i] = g->begin[i - 1] + deg[i - 1];
  } 
  
  printf("Adding edges\n");
  for (edge_t i = 0; i < M; i++) {
    node_t u = edges[i].first;
    node_t v = edges[i].second;

    edge_t pos = deg[u]--;
    assert(pos > 0);
    g->node_idx[g->begin[u] + pos - 1] = v;  // set end node of this edge
  }
  gettimeofday(&T2, NULL);
  printf("Manipulation time (ms) = %lf\n", ((T2.tv_sec) - (T1.tv_sec)) * 1000 + (T2.tv_usec - T1.tv_usec) * 0.001);

  printf("Storing binary\n");
  gettimeofday(&T1, NULL);
  g->store_binary(filename);
  gettimeofday(&T2, NULL);
  printf("storing time (ms) = %lf\n", ((T2.tv_sec) - (T1.tv_sec)) * 1000 + (T2.tv_usec - T1.tv_usec) * 0.001);
  printf("Done, freeing memory\n");
  delete g;
  // delete[] src;
  // delete[] dst;
  delete[] deg;
}

void write_matrix_market_file(char* filename, vector< pair<node_t, node_t> > edges, node_t N, edge_t M){
  FILE *ofp;
  if ((ofp = fopen(filename, "w")) == NULL) {
    printf("Error: output file is invalid\n");
    exit(1);
  }

  printf("Found %i Nodes, %i Edges\n", N, M);

  MM_typecode matcode;
  mm_initialize_typecode(&matcode);
  mm_set_matrix(&matcode);
  mm_set_coordinate(&matcode);
  mm_set_real(&matcode);

  mm_write_banner(ofp, matcode); 
  mm_write_mtx_crd_size(ofp, N, N, M);
  for(unsigned int i = 0; i<M; i++){
    fprintf(ofp, "%d %d 1.0\n", edges[i].first, edges[i].second);
  }

  fclose(ofp);
}

void print_analysis(vector< pair<node_t, node_t> > edges, node_t N, edge_t M){
  // allocate space for edges, and degree counts
  edge_t* deg = new edge_t[N];
  edge_t* ine = new edge_t[N];
  edge_t max_degree = 0; node_t max_node = 0;
  edge_t tail_verticies = 0;
  edge_t head_verticies = 0;
  edge_t disc_verticies = 0;
  memset(deg, 0, sizeof(edge_t) * N);
  memset(ine, 0, sizeof(edge_t) * N);

  gettimeofday(&T1, NULL);
  printf("Reading from vector to dense data.\n");
  //assign to our internal structures
  for(unsigned int i = 0; i<M; ++i) //iterate over edges - should use iterator :/
  {
    deg[edges[i].first]++;
    ine[edges[i].second]++;
    if(deg[edges[i].first] > max_degree){
    	max_degree = deg[edges[i].first];
    	max_node = edges[i].first;
    }
  }
  for(unsigned int i = 0;i<N; ++i){
  	if(deg[i] == 0 && ine[i] != 0){
  		tail_verticies++;
  	}
    if(ine[i] == 0 && deg[i] != 0){
      head_verticies++;
    }
    if(ine[i] == 0 && deg[i] == 0){
      disc_verticies++;
    }
  }
  printf("Max degree: %i from node %i\n", max_degree, max_node);
  printf("Tail verticies: %i\n", tail_verticies);
  printf("Head verticies: %i\n", head_verticies);
  printf("Disconnected verticies: %i\n", disc_verticies);
  float sparseSize = ((float)(N*max_degree))/((float)(4*1024*1024));
  printf("Unragged sparse size: %fGB\n", sparseSize);
  float denseSize = ((float)(N*N))/((float)(4*1024*1024));
  printf("Dense size: %fGB\n", denseSize);
  printf("\n\n");
  printf("Usable for GTX480 dense:  %s\n", (denseSize<1.53)?"yes":"no");
  printf("Usable for GTX480 sparse: %s\n", (sparseSize<1.53)?"yes":"no");
  printf("Usable for    K40 dense:  %s\n", (denseSize<11.5)?"yes":"no");
  printf("Usable for    K40 sparse: %s\n", (sparseSize<11.5)?"yes":"no");

  delete[] deg;
}

int main(int argc, char** argv){
  vector< pair<node_t, node_t> > edges;
  node_t vertex_count; 
  edge_t edge_count;
  bool undirected = false;
  bool analyse = false;
  char* g_in = NULL;
  char* gm_out = NULL;
  char* mm_out = NULL;

  for(int i = 1;i<argc;i++) //iterate over arguments - the first arg won't ever be one we care about.
  {
    if(strcmp(argv[i], "-u") == 0){ // parse an undirected file
      g_in = argv[i+1];
      undirected = true;
    }
    if(strcmp(argv[i], "-d") == 0){ // parse a directed file
      g_in = argv[i+1];
      undirected = false;
    }
    if (strcmp(argv[i], "-g") == 0){ //write to a green marl file
      gm_out = argv[i+1];
    }
    if (strcmp(argv[i], "-m") == 0){ //write to a green marl file
      mm_out = argv[i+1];
    }
    if (strcmp(argv[i], "-a") == 0){ //write to a green marl file
      analyse = true;
    }
  }
  if(g_in == NULL){
    printf("No input graph file specified. Failing.\n");
    printf("Specify a file using either -u (for undirected) or -d (for directed)\n");
    exit(1);
  }
  if((gm_out == NULL) and (mm_out == NULL) and (!analyse)){
    printf("No output graph file specified, and analysis not expected. Failing.\n");
    printf("Specify a file using either -g (for green-marl) or -m (for matrix-market)\n");
    exit(1);
  }
  printf("Input file %s\n", g_in);
  edges = parse_adjacency_file(g_in, undirected);
  vertex_count = max_node(edges)+1;
  edge_count = edges.size();

  if(gm_out != NULL){
    write_green_marl_file(gm_out, edges, vertex_count, edge_count);
  }
  if(mm_out != NULL){
    write_matrix_market_file(mm_out, edges, vertex_count, edge_count);
  }
  if(analyse){
    print_analysis(edges, vertex_count, edge_count);
  }

  return 0;
}
