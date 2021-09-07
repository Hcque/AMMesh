#include "Algorithm.h"
#include <algorithm>
#include <queue>
#include <utility>
#include <vector>
#include <cassert>

typedef std::pair<int, int> P;
static std::vector<int> dist, vis, fa;
static const int INF = 1e9;

static inline double Caldis(MVert* a, MVert* b){
	return std::sqrt( 
		(a->x()-b->x())* (a->x()-b->x())
		+(a->y()-b->y())* (a->y()-b->y())
		+(a->z()-b->z())* (a->z()-b->z())
	);
}

//��������lmk֮��ĵ��·��
void Dijkstra(PolyMesh& Mesh, int s_p, int e_p, std::vector<int>& path, int& tot)
{
	int n = Mesh.numVertices();
	dist.resize(n); vis.resize(n); fa.resize(n);
	for (int i=0;i<n;i++) {
		dist[i] = INF; vis[i] = 0; fa[i] = -1;
	}

	std::priority_queue<P> que; 
	que.push(std::make_pair(0, s_p) );
	dist[s_p] = 0;
	while (!que.empty()){
		P now = que.top(); que.pop();
		MVert* v = Mesh.vert(now.second);
		for (MVert* nxt: Mesh.vertAdjacentVertices(v)){
			int nxtId = nxt->index();
			double w = Caldis(v, nxt);
			if ( dist[now.second] + w < dist[nxtId]){
				dist[nxtId] = dist[now.second] + w;
				fa[nxtId] = now.second;
				que.push(std::make_pair(-w,nxtId));
			}
		}
	}
	int x = e_p;
	std::vector<int> revPath;
	while (x!= s_p) {
		assert(x >= 0);
		revPath.push_back(x);
		x = fa[x];
	}
	revPath.push_back(s_p);
	for (int i=0;i<revPath.size();i++){
		path.push_back(revPath[ revPath.size()-i-1]);
	}
	tot = dist[e_p];
}

void Dijkstra(PolyMesh& Mesh, int s_p, int e_p, std::vector<int>& path){
	int tmp = -1;
	Dijkstra(Mesh, s_p, e_p, path, tmp);
}
int Find(int x){
	if (fa[x] == x) return x;
	return fa[x] = Find(fa[x]);
}
bool Same(int x, int y){
	int u = Find(x); int v = Find(y);
	return (u==v);
}
void Union(int x, int y){
	int u = Find(x); int v = Find(y);
	if (u==v) return;
	fa[u] = v;
}
void Dijkstra_group(PolyMesh& Mesh, std::vector<int>& lmk, std::vector<std::vector<int>>& path)
{
	int n = lmk.size();
	std::vector<pathInfo> E;
	// foyld
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++){
			if (i==j) continue;
			int tot = 0;
			std::vector<int> _path;
			Dijkstra(Mesh, lmk[i], lmk[j], _path, tot);
			E.push_back(pathInfo(i+1, j+1, tot, _path)); // vert 1, ... n
		}
	}
	// mst
	fa.clear(); fa.resize(n);
	for (int i=1;i<=n;i++) fa[i] = i;

	std::sort(E.begin(), E.end());
	for (int i=0;i<E.size();i++){
		if ( Same(E[i].from, E[i].to)) continue;
		Union(E[i].from, E[i].to);
		E[i].selected = 1;
	}

	path.clear();
	for (int i=0;i<E.size();i++){
		if (E[i].selected){
			path.push_back(E[i].path);
			// set color  of chosen path
			for (int j=0;j<E[i].path.size();j++){
				Mesh.vert(E[i].path[j])->setColor(50,50,50);
			}
		}
	}
}