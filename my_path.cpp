#include <iostream>
#include <queue>
#include "Eigen/Dense"
#include <vector>
#include <chrono>
#include <stack>

using namespace std;
//time measurement
std::stack<clock_t> tictoc_stack;
void tic();
void toc();

int inf = 999999999; //no graph edge weight is more than this

priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> min_unv_costs; //holds the costs of all unvisited nodes that have been measured, and node # 

typedef pair<int, int> dest_weight; //used for adjacency list 

int lowest_unv_cost(const Eigen::MatrixXi &unv, const Eigen::MatrixXi &vis, int col_num);

vector<int> dijkstra (int src, int dest, vector<vector<dest_weight>>  &graph,
                      Eigen::MatrixXi &unv_ls, Eigen::MatrixXi &vis_ls)
{
    using namespace Eigen;

    //edge-case: src = dest
    if (src == dest)
        return {-1};

    vector<int> shortest_path {0};

    unv_ls.col(0).array() = inf;
    unv_ls.col(1).array() = -2; //-2 means unknown / no entry

    // -2 means no entry / unvisited
    vis_ls.col(0).array() = -2;
    vis_ls.col(1).array() = -2;
    
    //init src node to 0 cost and add to pqueue to be explored first
    unv_ls(src, 0) = 0;  
    unv_ls(src, 1) = -1; //prev of -1 means this node is source

    min_unv_costs.push(make_pair(0,src)); //pqueue insertion is slow so it should not replace unv_ls matrix

    //if we have visited the destination and added it to the vis_ls, stop searching
    bool vis_dest = false;
    while (!vis_dest)
    {  
        int visiting_node =  min_unv_costs.top().second;
        while (vis_ls(visiting_node,0) != -2)
        {
            min_unv_costs.pop();
            visiting_node =  min_unv_costs.top().second;
        }

        for (int i = 0; i < graph[visiting_node].size(); i++)
        {
            // cycle trough the adj list, at index of visiting node, which itself holds a vector of all adj nodes
            if (vis_ls(graph[visiting_node][i].first,0) == -2) // this node is not yet visited (dont visit nodes alr visited,
                                                               // as they alr have best shortest path)
            {
                // if edge between visiting_node & adjacent node + unvisited list cost of visiting node
                // is < cost of adjacent node's unvisited list cost, overwrite unvisited list (adjacent node)
                // with left side of comparison

                if (graph[visiting_node][i].second + unv_ls(visiting_node, 0)
                    < unv_ls(graph[visiting_node][i].first, 0))
                {
                    //overwrite cost, previous of adjacent node of visited node in unvisited list
                    unv_ls (graph[visiting_node][i].first, 0) = graph[visiting_node][i].second + unv_ls(visiting_node, 0);
                    unv_ls (graph[visiting_node][i].first, 1) = visiting_node;

                    //add to new cost to pqueue (cost, node)
                    min_unv_costs.push(make_pair(graph[visiting_node][i].second + unv_ls(visiting_node, 0),
                                                 graph[visiting_node][i].first));
                } //else do nothing (do not overwrite)
            }
        }
        // done visiting this node, now add it to visited list
        vis_ls(visiting_node,0) = unv_ls(visiting_node,0);
        vis_ls(visiting_node,1) = unv_ls(visiting_node,1);

        if (visiting_node == dest)
        {
            int prev = visiting_node;
            shortest_path[0] = prev;
            while(prev != src)
            {
                shortest_path.push_back(vis_ls(prev,1));
                prev = vis_ls(prev,1);
            }
            vis_dest = true;
        }
           
    }
    return shortest_path;
}

//returns row # of matrix with lowest value in col given
int lowest_unv_cost(const Eigen::MatrixXi &unv, const Eigen::MatrixXi &vis, int col_num)
{
    int lowest_index = -1;
    int lowest_found = inf;
    for (int i = 0; i < unv.rows(); i++)
    {
        if(vis(i,0) == -2) //only consider unvisited nodes
        {
            if (unv(i,col_num) < lowest_found)
            {
                lowest_index = i;
                lowest_found = unv(i,col_num);
            }
        }
    }
    return lowest_index;
}

int main (void)
{
    using namespace Eigen;

    int n = 5; //# of nodes
    //index of matrix is node (EX: node B is matrix.row(1))
    MatrixXi unv_ls(n, 3); //unvisited list, cols: | cost/g-score | previous |
    MatrixXi vis_ls(n, 2); //visited list, cols: | cost/g-score | previous |
    vector<vector<dest_weight>> graph(n); //vec of size n

    //cout << unv_ls.rows() << " " << unv_ls.cols() << endl;

    //------------------- Test Graph Construction (adj list) -----------------
    // Node A
    graph[0].push_back({1, 8});
    graph[0].push_back({2, 5});

    //B
    graph[1].push_back({0, 8});
    graph[1].push_back({2, 3});
    graph[1].push_back({3, 1});

    //C
    graph[2].push_back({0, 5});
    graph[2].push_back({1, 3});
    graph[2].push_back({3, 6});
    graph[2].push_back({4, 9});

    //D
    graph[3].push_back({1, 1});
    graph[3].push_back({2, 6});
    graph[3].push_back({4, 2});

    //E
    graph[4].push_back({2, 9});
    graph[4].push_back({3, 2});
 
    //--------------------------------------------------------------------
    // min_unv_costs.push(make_pair(2,20));
    // min_unv_costs.push(make_pair(8,50));
    // cout << min_unv_costs.top().first << "  " << min_unv_costs.top().second << endl;

    auto start = chrono::steady_clock::now();
    vector<int> shortest_path = dijkstra(0, 4, graph, unv_ls, vis_ls);
    auto end = chrono::steady_clock::now();

    for (int i = 0; i < shortest_path.size(); i++)
    {
        cout << shortest_path[i];
        if (i != shortest_path.size() - 1)
            cout << " --> ";
    }
    cout << endl;


    auto diff = end - start;
    cout << "Time: " << diff.count()/1000 << " micro sec" << endl;
    return 0;

}
// time measurement functions
void tic()
{
    tictoc_stack.push(clock());
}

void toc()
{
    cout.precision(10);
    std::cout << "Time elapsed: "
              << ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
              << std::endl;
    tictoc_stack.pop();
}