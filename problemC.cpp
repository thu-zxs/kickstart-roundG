#include <iostream>
#include <vector>
#include <fstream>
// #include <queue>
using namespace std;

const int INF = 100000000;

class Node
{
	public:
		int index;
		// int d;
		int value;
		Node* next;
		Node():next(NULL) {}
};


void Dijkstra( Node* &Adj, int src, vector<int> &shortest_dist, vector<int> &S, int frozen)
{
	int max_size = shortest_dist.size();
	shortest_dist[src] = 0; // S[src] = 1;
	for (int i = frozen; i < max_size; i++) {
		int tmp = INF;
		int u = -1;
		for (int j = 0; j < max_size; j++) {
			if (!S[j] && shortest_dist[j] < tmp) {
				tmp = shortest_dist[j];
				u = j;
			}
		}
		if (u == -1) break;
		S[u] = 1;
		Node* p = Adj[u].next;
		while(p != NULL) {
			if (shortest_dist[u] + p->value < shortest_dist[p->index]) {
				shortest_dist[p->index] = shortest_dist[u] + p->value;
			}
			p = p->next;
		}
	}
}

int main() {
	int T;
	ifstream fin("inc.txt");
	ofstream fout("outc.txt");
	fin >> T;
	for (int t = 0; t < T; t++) {
		int N, M, E, Sr, Sc, Tr, Tc;
		fin >> N >> M >> E >> Sr >> Sc >> Tr >> Tc;
		Sr--; Sc--; Tr--; Tc--;
		int **V = new int* [N];
		for (int row = 0; row < N; row++) {
			V[row] = new int [M];
			for (int col = 0; col < M; col++) {
				fin >> V[row][col];
				V[row][col] = -V[row][col];
			}
		}
		// Adjacent list construction
		Node *Adj = new Node[N*M];
		vector<int> S(N*M, 0); // whether a node is in set `s`.
		int frozen = 0;
		int src = Sr*M + Sc;
		int end = Tr*M + Tc;
		for (int row = 0; row < N; row++) {
			for (int col = 0; col < M; col++) {
				Node* vertex = &Adj[col+row*M];
				if (V[row][col] == 100000) {frozen += 1; S[col+row*M] = 1; continue;}
				// l, u, r, d
				if (col-1 >= 0 && V[row][col-1] < 100000) {
					vertex->next = new Node();
					vertex->next->index = (col-1)+row*M;
					vertex->next->value = V[row][col-1];
					// vertex->next->d = INF;
					vertex = vertex->next;
				}
				if (row-1 >= 0 && V[row-1][col] < 100000) {
					vertex->next = new Node();
					vertex->next->index = (col)+(row-1)*M;
					vertex->next->value = V[row-1][col];
					// vertex->next->d = INF;
					vertex = vertex->next;
				}
				if (col+1 < M && V[row][col+1] < 100000) {
					vertex->next = new Node();
					vertex->next->index = (col+1)+(row)*M;
					vertex->next->value = V[row][col+1];
					// vertex->next->d = INF;
					vertex = vertex->next;
				}
				if (row+1 < N && V[row+1][col] < 100000) {
					vertex->next = new Node();
					vertex->next->index = (col)+(row+1)*M;
					vertex->next->value = V[row+1][col];
					// vertex->next->d = INF;
					vertex = vertex->next;
				}
			}
		}
		// priority_queue<int, vector<int>, greater<int> > shortest_dist(N*M, INF);
		vector<int> shortest_dist(N*M, INF);
		Dijkstra(Adj, src, shortest_dist, S, frozen);
		
		int ans = E >= shortest_dist[end] ? E-shortest_dist[end]:-1;
		fout << "Case #" << t+1 << ": " << ans << endl;
	}
	return 0;
}
