#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/random.hpp"
#include "boost/graph/make_connected.hpp"
#include "boost/random/mersenne_twister.hpp"
#include <boost/graph/grid_graph.hpp>

//For random number generator
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_int_distribution.hpp"
#include "boost/random/variate_generator.hpp"

#include <fstream>
#include <stack>
#include <functional>
#include <queue>
#include <vector>
#include <iostream>
#include <random>
#include <ctime>

using namespace boost;
using namespace std;

struct EdgeProperties 
{
    int weight; //edge weight property
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, EdgeProperties> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;

//Iterators

typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;
typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIterator;
typedef boost::property_map<Graph, int EdgeProperties::*>::type WeightMap;
typedef std::pair<int, int>E;

int columns; //global variable for grid columns

//Create random numbers from [min, max]
int randomRange(int min, int max)
{
    struct timespec tm;
    clock_gettime(CLOCK_REALTIME, &tm);
    double time = tm.tv_nsec;
    
    random::mt19937 rng(time);
    random::uniform_int_distribution<> six(min,max);
    
    return six(rng);
}

//Greate a rows x columns grid graph with edge weights [cost_min, cost_max]
Graph createGridGraph(int rows, int columns, int cost_min, int cost_max)
{
    int n_vertices = rows * columns;
    EdgeIterator ei, ei_end;
    Graph G(n_vertices);
   
    for (unsigned int i = 0; i < n_vertices; ++i) {
         
         if (i%n_vertices == n_vertices-1) 
            continue;
        
        if  ( (i + 1) % columns != 0 && i != n_vertices-1) {
            add_edge(i, i+1, G);
            add_edge(i+1, i, G);
        }

        if (i < n_vertices - columns) {
            add_edge(i, i+columns, G);
            add_edge(i+columns, i, G);        
        }
    }      
    // Converting to EdgeList graph
    vector<pair<Vertex, Vertex> > pairlist;
    for (tie(ei, ei_end) = edges(G); ei!=ei_end; ++ei)
        pairlist.emplace_back(source(*ei, G), target(*ei, G));

    Graph g(pairlist.begin(), pairlist.end(), num_vertices(G));

    //Add the weights to edges
    int random_weight;
    for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
         random_weight = randomRange(cost_min, cost_max);
        g[*ei].weight = random_weight;
         g[*++ei].weight = random_weight;
    }
    return g;
}

//Print grid graph to the file
void printGraphVizToFile(Graph g, string filename)
{
    ofstream gout;
    gout.open(filename);

    gout << "digraph G" << endl;
    gout << "{" << endl;

    VertexIterator v, v_end;
    for(tie(v, v_end) = vertices(g); v != v_end; ++v) {
        gout << "   " << *v << ";" << endl;
    }

    EdgeIterator e, e_end;
    for(tie(e, e_end) = edges(g); e != e_end; ++e) {
        gout << "   " << source(*e, g) << " -> " << target(*e, g);
        gout << " [label=" << g[*e].weight << "];" << endl;
    }
    gout << "}" << endl;
}

void setColumns(int c){
    columns = c;
}

int getColumns() {
    return columns;
}

bool dijkstra_shortest(const Graph& G, const int& s, const WeightMap& weight, vector<long>& dist, bool show_info) { 
    typedef pair<Vertex, long> pairQueue;    
    int N = num_vertices(G);
    priority_queue< pairQueue, vector <pairQueue>, greater<pairQueue> > pq; 
  
    const Vertex NIL =  graph_traits<Graph>::null_vertex(); // end marker

    dist[s] = 0;
    pq.push(make_pair(s, dist[s]));

    Vertex v, u;
    OutEdgeIterator oei, oei_end;
    
    while (!pq.empty())  //until priority queue is not empty
    {   
        u = pq.top().first;
        pq.pop(); 
  
        // Get all adjacent of u.  
        for (tie(oei, oei_end) = out_edges(u, G); oei != oei_end; ++oei) {
            
            int v = target(*oei, G); //proskimeni tou u 
            int weight_adj = weight[*oei]; 
            int alt = dist[u] + weight_adj;
            
            if (dist[v] > alt)
            { 
                dist[v] = alt;
                pq.push(make_pair(v, dist[v]));
            }
        }
    }
    return true;
}

bool dijkstra_sp(const Graph& G, const int& s, const int& t, const WeightMap& weight, vector<long>& dist, bool show_info) { 
    typedef pair<Vertex, long> pairQueue;    
    
    int N = num_vertices(G);
    int nodes_visited = 0;

    priority_queue< pairQueue, vector <pairQueue>, greater<pairQueue> > pq; 
    const Vertex NIL =  graph_traits<Graph>::null_vertex(); // end marker

    dist[s] = 0;
    pq.push(make_pair(s, dist[s]));

    Vertex v, u;
    OutEdgeIterator oei, oei_end;
    
    while (!pq.empty()) 
    {   
        u = pq.top().first; 
        pq.pop();

        if (u == t) { //until t will pop from priority queue then end
            break;
        }

        nodes_visited++;
        // Get all adjacents of u.  
        for (tie(oei, oei_end) = out_edges(u, G); oei != oei_end; ++oei) {
            
            int v = target(*oei, G); //proskimeni tou u 
            int weight_adj = weight[*oei]; 
            int alt = dist[u] + weight_adj;
            if (dist[v] > alt)
            { 
                dist[v] = alt;
                pq.push(make_pair(v, dist[v]));
                
            }
        }
    }
    
    if (show_info) {
        cout << endl << "************** Dijkstra-SP Results **************" << endl;
        cout << "Vertex   Distance from Source\n"; 
        for (int i = 0; i < N; ++i) {
            cout <<  i << "\t" << dist[i] << endl; 
        }
    }
    cout << "Dijkstra-SP nodes visited: " << nodes_visited << endl;
    return true;
}

bool a_star(const Graph& G, const int& s, const int& t, const WeightMap& weight, vector<long>& dist, bool show_info) { 
    typedef pair<Vertex, long> pairQueue;

    int N = num_vertices(G);
    int nodes_visited = 0;

    priority_queue< pairQueue, vector <pairQueue>, greater<pairQueue> > pq; 
  
    const Vertex NIL =  graph_traits<Graph>::null_vertex(); // end marker

    dist[s] = 0;
    int x_s  = s / getColumns();
    int y_s = (s % getColumns());
    pq.push(make_pair(s, dist[s]));
    Vertex v, u;
    OutEdgeIterator oei, oei_end;

    std::vector<Vertex> shortest_path;
    int h_t = 0;
    int x_t = t / getColumns();
    int y_t = v / getColumns();
    int h_s = sqrt( (x_s - x_t)^2 + (y_s - y_t)^2 );
    while (!pq.empty()) 
    {   
        u = pq.top().first; 
        shortest_path.push_back(u);
        pq.pop(); 
  
        // Get all adjacent of u.  
        for (tie(oei, oei_end) = out_edges(u, G); oei != oei_end; ++oei) {
            
            int v = target(*oei, G); //proskimeni tou u  
           int weight_adj = weight[*oei]; 
            
            int x_i = u / getColumns();
            int x_j = v / getColumns();
            int y_i = (u % getColumns());
            int y_j = (v % getColumns());
            int h_i = sqrt( (x_i - x_t)^2 + (y_i - y_t)^2 );
            int h_j = sqrt( (x_j - x_t)^2 + (y_j - y_t)^2 );
           
            int new_weight = weight_adj + h_j - h_i; 
            int alt = dist[u] + new_weight;
            nodes_visited++;
            if (dist[v] > alt)
            { 
                dist[v] = alt;
                pq.push(make_pair(v, dist[v]));
            }
        }

        if (u == t) {
            dist[u] = dist[u] - (h_t - h_s);
            break;
        }
    }

    if (show_info) {
        cout << endl << "************** A* Results **************" << endl;
        cout << endl << "Vertex   Distance from Source\n"; 
        for (int i = 0; i < N; ++i) {
            cout <<  i << "\t" << dist[i] << endl;  
        }
        //cout << endl;
        //for (auto i = shortest_path.begin(); i != shortest_path.end(); ++i) 
          //  cout << *i << " "; 
        //cout << endl;
    }

    cout << "A* nodes visited: " << nodes_visited << endl;
    return true;
}

int main(int argc, char *argv[])
{
    bool show_info = false;

     if (argc < 5) {
        cout << "Specify number of rows, columns, min and max weights" << endl;
        return 0;
    }
    else if (argc == 5)
         show_info = false;
    else if (argc == 6 && std::string(argv[5])=="info")
        show_info = true;
    else
        return 0;

    const int rows = atoi(argv[1]);
    const int columns = atoi(argv[2]);
    setColumns(columns);
    const int min_weight = atoi(argv[3]);
    const int max_weight = atoi(argv[4]);

    const int n_vertices = rows * getColumns();
    const int infinite = (numeric_limits < int >::max)();

    Graph g = createGridGraph(rows, getColumns(), min_weight, max_weight);
    EdgeIterator ei, ei_end;
    int i = 0;
    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap mIndexMap = get(vertex_index, g);
    WeightMap weight = get(&EdgeProperties::weight, g);
    pair<VertexIterator, VertexIterator> vip;

    if (show_info) {
        cout << "Komboi:\n"; 
        for(vip = vertices(g); vip.first != vip.second; ++vip.first){
            cout << mIndexMap[*vip.first] << " ";
        }
        
        cout << endl << endl << "Akmes:\n";
        for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            cout << "(" << mIndexMap[source(*ei, g)] << ", "
                    << mIndexMap[target(*ei, g)] << ") " << g[*ei].weight << ", " ;
        }
        cout << endl;
    }

    //creating s
    int random_columns_s[rows];
    int j = 0;

    for (unsigned int i = 0; i < n_vertices; i++) {
        if (i % getColumns() == 0) {
            random_columns_s[j] = i;
            j++;
            if (j==rows)
            break;
        } 
    }

    //creating t
    int random_columns_t[rows];
    j = 0;

    for (i = 0; i < n_vertices; i++) {
        if ( (i+2) % getColumns() == 0 ) {
            random_columns_t[j] = i;
            j++;
        }
    }

    int random_column_index = randomRange(0, rows-1);
    int s = random_columns_s[random_column_index];
    int random_column_idx = randomRange(0, rows-1);
    int t = random_columns_t[random_column_idx];
    
    if (t < s) {
        int temp = s;
        s = t;
        t = temp;
    }

    cout << "s: " << s << endl;
    cout << "t: " << t << endl;
    
    vector<long> dist(n_vertices, infinite);
    printGraphVizToFile(g, "graph.dot");
    
    int exec_times = 1;

    clock_t begin = clock();

    for (int i = 0; i < exec_times; ++i)
        bool result_dijkstra = dijkstra_shortest(g, s, weight, dist, show_info);
    clock_t end = clock();
    double elapsed_secs_dijkstra = double(end - begin) / CLOCKS_PER_SEC / exec_times;

    begin = clock();
    for (int i = 0; i < exec_times; ++i)
        bool result_dijkstra_sp = dijkstra_sp(g, s, t, weight, dist, show_info);
    end = clock();
    double elapsed_secs_dijkstra_sp = double(end - begin) / CLOCKS_PER_SEC / exec_times;

    begin = clock();
    for (int i = 0; i < exec_times; ++i)
        bool result_a_star = a_star(g, s, t, weight, dist, show_info);
    end = clock();
    double elapsed_secs_a_star = double(end - begin) / CLOCKS_PER_SEC / exec_times;

    cout << endl << "Benchmark results: " << endl;
    cout << "Time elapsed Dijkstra: " << elapsed_secs_dijkstra << " seconds" << endl;
    cout << "Time elapsed Dijkstra-SP: " << elapsed_secs_dijkstra_sp << " seconds" << endl;
    cout << "Time elapsed A*: " << elapsed_secs_a_star << " seconds" << endl;
    
    return 0;
}
