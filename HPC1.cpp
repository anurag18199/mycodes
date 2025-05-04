#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <omp.h>

using namespace std;

class Graph {
    int V;
    vector<vector<int>> adj;

public:
    Graph(int V) {
        this->V = V;
        adj.resize(V);
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u); // Undirected graph
    }

    void parallelBFS(int start) {
        vector<bool> visited(V, false);
        queue<int> q;

        visited[start] = true;
        q.push(start);

        cout << "Parallel BFS: ";

        while (!q.empty()) {
            int level_size = q.size();

            // Process current level
            #pragma omp parallel for
            for (int i = 0; i < level_size; i++) {
                int node;
                
                #pragma omp critical
                {
                    node = q.front();
                    q.pop();
                }

                cout << node << " ";

                // Visit all neighbors
                for (int neighbor : adj[node]) {
                    if (!visited[neighbor]) {
                        #pragma omp critical
                        {
                            if (!visited[neighbor]) {
                                visited[neighbor] = true;
                                q.push(neighbor);
                            }
                        }
                    }
                }
            }
        }

        cout << endl;
    }

    void parallelDFSUtil(int node, vector<bool>& visited) {
        visited[node] = true;
        cout << node << " ";

        #pragma omp parallel for
        for (int i = 0; i < adj[node].size(); i++) {
            int neighbor = adj[node][i];
            if (!visited[neighbor]) {
                #pragma omp task
                parallelDFSUtil(neighbor, visited);
            }
        }

        #pragma omp taskwait
    }

    void parallelDFS(int start) {
        vector<bool> visited(V, false);
        cout << "Parallel DFS: ";

        #pragma omp parallel
        {
            #pragma omp single
            parallelDFSUtil(start, visited);
        }

        cout << endl;
    }
};

int main() {
    Graph g(7);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 3);
    g.addEdge(1, 4);
    g.addEdge(2, 5);
    g.addEdge(2, 6);
    double t1, t2;
    g.parallelBFS(0);
    // t1 = omp_get_wtime();
    // g.parallelBFS(0);
    // t2 = omp_get_wtime();
    // cout << "Sequential Merge Sort:  " << (t2 - t1) << " s\n";
    
    // t1 = omp_get_wtime();
    // g.parallelDFS(0);
    // t2 = omp_get_wtime();
    // cout << "Sequential Merge Sort:  " << (t2 - t1) << " s\n";
    

    return 0;
}
