// mesh simplication 
// ===============================================

// caution summary:

// set <opertoar
// insert to set , const &
// local for faces
// idx needs repad when writing
// costomized heap, Hash

// int ? double !!!
// points. not index for counting !
// == / = !!
// local g variable


#include "PolyMesh/IOManager.h"
#include "Eigen/Dense"
#include <queue>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <time.h>
#include <chrono>

#define SHRAP_COST 1
#define DFS_SEARCH 1

using namespace acamcad;
using namespace polymesh;


int iters;
double t;
double simp_percentage ;
int initNumV;
// std::unordered_map<ValidPair, int, ValidPair::HashFunction> pairDup;
std::vector<int> toDel;

// std::unordered_map<int , Eigen::Matrix4d> findQ;
// Heap heap;


class Mesh
{
private:

	struct Face{
		mutable int dim[3];
		Face(int _i, int _j, int _k){
			dim[0] = _i;
			dim[1] = _j;
			dim[2] = _k;
			sort();
		}
		void sort(){
			int min_node = std::min( dim[0], std::min(dim[1], dim[2]));
			if (min_node == dim[0]) return;
			else if (min_node == dim[1]) {
				// 012 ->120
				std::swap(dim[0], dim[1]); std::swap(dim[1], dim[2]);
			}
			else if (min_node == dim[2])
			{
				// 012 -> 201
				std::swap(dim[0], dim[1]); std::swap(dim[0], dim[2]);
			}
			return;
		}
		// MVector3 getnormal()
		// {


		// }


	
		bool contain(int vi) const { 
			return std::find(dim,dim+3, vi) != dim+3;
		}

		int find_diff(int v1) const
		{
			if (dim[0] != v1) return dim[0];
			else if (dim[1] != v1) return dim[1];
			else if (dim[2] != v1) return dim[2];
			else assert(0);
		}
		int find_diff(int v1, int v2) const
		{
			std::set<int> copy(std::begin(dim), std::end(dim));
			copy.erase(v1);
			copy.erase(v2);
			return *copy.begin();
		}

		void replace(int v1, int v2)
		{
			if (dim[0] == v1) dim[0] = v2;
			else if (dim[1] == v1) dim[1] = v2;
			else if (dim[2] == v1) dim[2] = v2;
			else assert(0);
			return;
		}

		bool reset(int v2, int v1)const  // caution!
		{
			auto it = std::find(dim,dim+3, v2);
			if (it == dim+3) return false;
			dim[it - dim] = v1;
			return true;
		}

		bool find_other2(int v1,int& o2,int &o3) const
		{
			if (dim[0] == v1) o2 = dim[1], o3 = dim[2];
			else if (dim[1] == v1) o2 = dim[0], o3 = dim[2];
			else if (dim[2] == v1) o2 = dim[0], o3 = dim[1]; 
			// else assert(0);
			else return false;
			return true;
		}

		friend std::ostream& operator<<(std::ostream &out, const Face& f)
		{
			out << f.dim[0] << "|" << f.dim[1] << "|" << f.dim[2] << "\n";
			return out;
		}

		friend bool operator<(const Face &a, const Face& b)  
		{
			return (a.dim[0] != b.dim[0]) ? (a.dim[0] < b.dim[0]) :
			( (a.dim[1] != b.dim[1]) ? (a.dim[1] < b.dim[1]) :
			(a.dim[2] < b.dim[2]) );
		}
		bool operator==(const Face& other) const {
			return dim[0] == other.dim[0] &&
			 dim[1] == other.dim[1] &&
			 dim[2] == other.dim[2] ;
		} 

		struct HashFunction
		{
			std::size_t operator()(const Face& k) const
			{
			using std::size_t;
			using std::hash;

			// Compute individual hash values for first,
			// second and third and combine them using XOR
			// and bit shifting:

			return (hash<int>()(k.dim[0]*(13+k.dim[1]) ) ^ (hash<int>()(( k.dim[2] + 7) * k.dim[0]) << 1) );
			}
		};
	};

	struct Vertice
	{
		MPoint3 position;
		std:: set<Face> fas;
		Eigen::Matrix4d Q;
		int index;
		Vertice(){}
		Vertice(double i,double j, double k, int idx){
			position[0] = i;
			position[1] = j;
			position[2] = k;
			index = idx;
		}

		void replace(double i, double j, double k, Eigen::Matrix4d _Q)
		{
			position[0] = i;
			position[1] = j;
			position[2] = k;
			Q = _Q;

		}
		int getFaceNum()
		{
			return fas.size();
		}

		bool count_face(const Face& ff)
		{
			return fas.count(ff) > 0;
		}

		void add_face(const Face& f){ fas.insert(f); }
		void remove_face(const Face& f){ fas.erase(f); }

		double cal_dist(const Vertice& other)
		{
			return position.distance(other.position);
		}


		bool reset(int vnow, int vreplace)
		{
			bool ans = false;
			for (auto &ff: fas)
				ans |= ff.reset(vnow, vreplace);
			return ans;
		}

		std::set<int> oneRingNeibours()
		{
			std::set<int> ans;
			for (auto &ff: fas) 
				for (int x: ff.dim) if (x != index) ans.insert(x);
			return ans;
		}
	};

	struct ValidPair
	{
		// MEdge* me;
		int vi, vj, time;
		double cost;
		Eigen::Vector4d vhat;
		ValidPair() {}
		ValidPair(int _i, int _j): vi(_i), vj(_j){
			sort();
		}
		ValidPair(int _i, int _j, double _c, Eigen::Vector4d _hat): vi(_i), vj(_j),
							cost(_c) , vhat(_hat)  {
			sort();
		}
		bool operator<(const ValidPair& other) const { return cost > other.cost; } 
		bool operator==(const ValidPair& other) const { return vi == other.vi && vj == other.vj;  } 
		void sort()
		{
			if (vi > vj) std::swap(vi,vj);
		}

		unsigned long long index_hash() const {
			return ((unsigned long long)(vi) << 32) | (unsigned long long )(vj);
		}


	// // https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key

	//   struct HashFunction
	//   {
	//     std::size_t operator()(const ValidPair& k) const
	//     {
	//       using std::size_t;
	//       using std::hash;

	//       // Compute individual hash values for first,
	//       // second and third and combine them using XOR
	//       // and bit shifting:

	//       return (hash<int>()(k.vi) ^ (hash<int>()(k.vj) << 1) );
	//     }
	//   };

	};


	struct Heap
	{
		int tstamp;
		std::unordered_map<unsigned long long, int > times;
		std::priority_queue<ValidPair> que;
		std:: vector< std::set<int>> in_que_pairs;

		Heap() { tstamp = 0; }
		inline void resize_in_queue_pairs(int n)
		{
			in_que_pairs.resize(n);
		}
		void clear()
		{
			tstamp = 0; times.clear();
			while (que.size()) que.pop();
		}
		inline bool count(unsigned long long index)
		{
			return (times.count(index)) ? (times[index] > 0) : false;
		}

		inline int size() { return que.size(); }

		void refresh()
		{
			while (que.size())
			{
				auto pair = que.top();
				if (pair.time == times[pair.index_hash()]) break;
				que.pop();
			}
		}

		void pop()
		{
			refresh();
			que.pop();
		}

		ValidPair top()
		{
			refresh();
			return que.top();
		}

		void push(ValidPair& pair)
		{
			++tstamp ;
			pair.time = tstamp; times[pair.index_hash()] = tstamp;
			que.push(pair);
			in_que_pairs[pair.vi].insert(pair.vj);
			in_que_pairs[pair.vj].insert(pair.vi);
		}
		void del(const ValidPair& pair)
		{
			times[pair.index_hash()] = -1;
			in_que_pairs[pair.vi].erase(pair.vj);
			in_que_pairs[pair.vj].erase(pair.vi);
		}
	};

// class Heap : public std::priority_queue<ValidPair>
// {
// public:
// 	bool remove_related(const int v, int v2)
// 	{
// 		bool flag = false;
// 		for ( auto it = this->c.begin(); it != this->c.end(); ++it)
// 		{
// 			if ( (*it).vi == v || (*it).vj == v || (*it).vi == v2 || (*it).vj == v2){
// 				this->c.erase(it);
// 				flag = true;
// 			}	
// 		}
	
// 		std::make_heap(this->c.begin(), this->c.end(), this->comp);
// 		return flag;
// 	}
// };

public:
	Mesh() { 
		tot = tot_face = 0; vertices.clear(); 
		mean_edge_len = 0.0;
	}
	~Mesh() {}
	int tot, tot_face;
	std::vector<Vertice> vertices;
	std::vector<std::set<int>> edges;
	Heap heap;
	// for write faces
	std::vector<int> mapidx;
	std::unordered_set<Face, Face::HashFunction> fascnt;

	std::vector<int> visit;
	double mean_edge_len;

	int numFaces() ;
	int numVertices() const ;
	void loadMesh(std::string, Mesh* mesh);
	void writeMesh(std::string, Mesh* mesh);
	void simp(); // main procedure
	void cal_Q();
	void select_pairs();
	void contract(int v1, int v2, Eigen::Vector4d);
	void contract_v2(int v1, int v2, Eigen::Vector4d vhat);

	void add_kp(const Face& ff, Eigen::Matrix4d& Q);
	void addFace(int i,  Face);
	void addVertice(double i, double j, double k, int index);
	MVector3 getnormal(int v1, int v2, int v3);
	double cal_cost(int v1, int v2, Eigen::Vector4d& newx );
	void add_pair(int v0, int v1);
	std::set<int> oneRingNeibours(int vi);
	void dfs_search(int now, int);

};

int Mesh::numFaces() { 
	int ans = 0;
	for (int i = 0; i < vertices.size(); i ++ ) if (!toDel[i])
		ans += vertices[i].fas.size();
	return ans;
 }

int Mesh::numVertices() const { return vertices.size(); }
void Mesh::addFace(int i,  Face face)
{
	vertices[i].fas.insert(face);
}
void Mesh::addVertice(double i, double j, double k, int index)
{
	vertices.push_back(Vertice(i,j,k,index));
}
MVector3 Mesh::getnormal(int v1, int v2, int v3) 
{
	MVector3 ba = vertices[v2].position - vertices[v1].position;
	MVector3 ca = vertices[v3].position - vertices[v1].position;
	// for (int k = 0 ; k< 3; k ++ ) 
	// 	std::cout << vertices[v2].position[k] <<" ";
	// std::cout << "\n";
	ba.normalize(); ca.normalize();
	return ba.cross(ca);
}


std::set<int> Mesh::oneRingNeibours(int vi)
{
	std::set<int> ans;
	for (auto &ff: vertices[vi].fas) 
		for (int x: ff.dim) if (x != vi) ans.insert(x);
	return ans;
}


void Mesh::contract_v2(int v1, int v2, Eigen::Vector4d vhat)
	{
		auto faceset = vertices[v2].fas; // f1,f2, 5,6,7
		for (auto &ff : faceset) {
			if (ff.contain(v1)) { // f1,f2
				int vother = ff.find_diff(v1,v2);
				vertices[vother].remove_face(ff);
				vertices[v2].remove_face(ff);
				vertices[v1].remove_face(ff);
				tot_face--;
			}
			else{  // 5,6,7
			assert(ff.contain(v2));

				int o2,o3;
				ff.find_other2(v2,o2,o3);
				assert(o2 != v1);
				assert(o3 != v1);
				vertices[o2].remove_face(ff);
				vertices[o3].remove_face(ff);

				Face newf(v1,o2,o3);
				if (! vertices[v1].count_face(newf)){
				vertices[v1].add_face(newf);
				vertices[o2].add_face(newf);
				vertices[o3].add_face(newf);
				}
				else {
					--tot_face;
				}
			}
		}
		// std::cerr << "constraction" << v2 << "\n";
		// vertices
		toDel[v2] = 1; tot--;
		// connect neigbours  to v1
		// for (auto &sf : vertices[v1].fas) for (int j : sf.dim ) if (j != v1  && j < vertices.size() && !toDel[j])
		// 	vertices[j].reset(v2,v1);

		// update v1
		vertices[v1].replace(vhat[0], vhat[1], vhat[2], vertices[v2].Q);

		// update heap
		std::set<int > in_que_pairs1 = heap.in_que_pairs[v1];
		std::set<int > in_que_pairs2 = heap.in_que_pairs[v2];
		for (auto v: in_que_pairs2)
		{
			heap.del(ValidPair(v,v2));
			if ( !in_que_pairs1.count(v) && v != v1) add_pair(v,v1);
		}

		for (auto v: in_que_pairs1) if (v != v2) add_pair(v,v1);
	}


void Mesh::loadMesh(std::string mesh_path, Mesh* mesh)
{
	double a, b, c; int fa,fb,fc;
	char ch;

	std::cout << "load..." << "\n";
	std::ifstream in(mesh_path);
	char buf[256];

	while (in.getline(buf, 256)){

		std::stringstream ss(buf);
		int cnt = 0;
		ss >> ch;
		switch (ch)
		{
			case '#': continue; break;
			case 'v':
				ss >> a >> b >> c;
				// std::cout << a << b << c << "\n";
				mesh->addVertice(a,b,c, cnt ++ );
				tot++;
				break;
			case 'f':
				ss >> fa >> fb >> fc;
				fa--; fb--; fc--;
				Face face(fa,fb,fc);
				mesh->addFace(fa, face);
				mesh->addFace(fb, face);
				mesh->addFace(fc, face);
				break;
		}
	}
}



int post_v_cnt= 0;
int post_f_cnt = 0;

void Mesh::writeMesh(std::string mesh_path, Mesh* mesh)
{
	std::cout << "write..." << "\n";
	int cnt = 0;
	std::ofstream out(mesh_path);

	std::cerr << "============" << vertices.size() << "\n";
	for (int i=0;i<vertices.size(); i++){
		if (toDel[i]) {
			mapidx.push_back(-2);
			continue;
		}
		mapidx.push_back(cnt++);
		out << 'v' << " " <<
			vertices[i].position[0] << " " <<
			vertices[i].position[1] << " " <<
			vertices[i].position[2] << std::endl;
		post_v_cnt ++ ;
	}
	for (int i=0;i<vertices.size(); i++){
		if (toDel[i]) continue;

		for ( auto &ff : vertices[i].fas){
			if (fascnt.find(ff) != fascnt.end()) continue;
			else{
			fascnt.insert(ff);

			out << 'f' << " " <<
			mapidx[ff.dim[0]]+1 << " " <<
			mapidx[ff.dim[1]]+1 << " " <<
			mapidx[ff.dim[2]]+1 << std::endl;
			post_f_cnt ++ ; 
			}
		}
	}
}

void Mesh::add_kp(const Face& ff, Eigen::Matrix4d& Q)
{
	MVector3 normal  = getnormal(ff.dim[0], ff.dim[1], ff.dim[2]);
	//  for (int k = 0 ; k< 3; k ++ ) 
    //         std::cout << normal(k) ;
    //     std::cout << "\n";
	assert( normal.dot(normal) - 1.0 < 0.001);
	auto tmp = vertices[ff.dim[0]].position;

	Eigen::Vector4d X(tmp[0], tmp[1], tmp[2] , 1);
	double d = -1.0f * ( normal[0]*X[0] + normal[1]*X[1] + normal[2]*X[2] );
	Eigen::Matrix<double, 4,1> p(normal[0], normal[1], normal[2], d);
	Eigen::Matrix4d Qi = p * p.transpose();
	Q += Qi;
        // for (int k = 0 ; k< 4; k ++ ) for (int j = 0; j < 4; j ++ )
        //     std::cout << Qi(k,j) << " ";
        // std::cout << "\n";
}

void Mesh::cal_Q()
{
	std::cout << "cal Q, stored in map" << std::endl;
	
	for (int i=0;i< vertices.size(); i++) {
		int vid = i;
		vertices[i].Q.setZero();
		for (auto& ff: vertices[i].fas)
			add_kp(ff, vertices[i].Q);
        // for (int k = 0 ; k< 4; k ++ ) for (int j = 0; j < 4; j ++ )
        //     std::cout << vertices[i].Q(k,j) ;
        // std::cout << "\n";

		// std::cout << vid << "|" << Q.data << "\n";
	}
}


double Mesh::cal_cost(int v1, int v2, Eigen::Vector4d& newx )
{
	MPoint3 X = vertices[v1].position;
	MPoint3 Y = vertices[v2].position;


	auto Qnew =  (vertices[v1].Q + vertices[v2].Q);

	double det = Qnew(0,0) * (Qnew(1,1)*Qnew(2,2)-Qnew(1,2)*Qnew(1,2)) - \
					Qnew(0,1) * (Qnew(0,1)*Qnew(2,2) - Qnew(1,2)*Qnew(2,0)) +\
					Qnew(0,2) * (Qnew(2,1)*Qnew(1,0) - Qnew(1,1)*Qnew(2,0)) ;
	if (fabs(det) < 1e-12){

	Eigen::Vector4d X4(X[0],X[1], X[2], 1.0f);
	Eigen::Vector4d Y4(Y[0],Y[1], Y[2], 1.0f);
	Eigen::Vector4d xhat = X4;
	newx = xhat;
	}
	else {
		Eigen::Matrix4d Qnewtmp = Qnew;
		Qnewtmp(3,0)  = Qnewtmp(3,1) = Qnewtmp(3,2) = 0; Qnewtmp(3,3) = 1;
		Eigen::Vector4d b(0,0,0,1);
		auto x = Qnewtmp.colPivHouseholderQr().solve(b);
		newx = x;
	}

	double ans =  ( newx.transpose() * Qnew * newx );

	// add penality for sharp details
#ifdef SHRAP_COST
	std:: vector<MVector3> normals;
	for (auto &ff: vertices[v1].fas){
		if (ff.contain(v2)){
			 normals.push_back( getnormal(ff.dim[0],ff.dim[1],ff.dim[2]) );
		}
	}
	if (normals.size() == 2){
		double penalty = normals[0].dot(normals[1]);
		// std::cerr << "penalty:" <<penalty << "\n";
		ans / std::min(1.1, 1.1 + penalty);
		// ans *std::exp( std::min(1.0, 1.1 + penalty) ) ;
	}

#endif
	return ans;
}

void Mesh::select_pairs()
{
	std::cout << "collect valid pairs, push in to heap" << std::endl;

	mean_edge_len = 0; int edge_cnt = 0;
	edges.resize( vertices.size() );
	heap.resize_in_queue_pairs( vertices.size() );

	for (int i = 0; i <vertices.size(); i ++ ) 
		for (auto &ff: vertices[i].fas)
			for (int &j : ff.dim ) 
				if (j != i) {
					edges[i].insert(j);
					mean_edge_len += vertices[i].cal_dist(vertices[j]);
					edge_cnt ++;
				}

	mean_edge_len /= edge_cnt;
	std:: cerr << "mean edge len: " << mean_edge_len << "\n";
	for (int i = 0; i < vertices.size(); i ++ ){
	// std::cerr << i << "\n";
		for (auto &j : edges[i]) if (i < j)
			add_pair(i, j);
	}

#ifdef DFS_SEARCH
	visit.resize(initNumV);
	for (int i = 0; i < vertices.size(); i ++ )
	{
		std::fill(visit.begin(), visit.end(),0);
		// std::cerr << "ds: " << i << "\n";

		dfs_search(i,i);
	}
#endif
	std::cout <<  "heap size: " << heap.size() << "\n";
}

void Mesh::dfs_search(int now, int origin)
{
	if (vertices[now].cal_dist(vertices[origin]) > mean_edge_len* t) return;

	if (now != origin) add_pair(origin, now);
	visit[now] = 1;
	for (auto i: vertices[now].oneRingNeibours())
	{
		// std::cerr << i << now <<  "|" << origin << "\n";
		if (!visit[i])
			dfs_search(i, origin);
	}
}

void Mesh::simp()
{
	// cal Q, stored in map
	cal_Q();
	// collect valid pairs, push in ti heap;
	select_pairs();

	// iters
	// while ( heap.size() && ( mesh->numVertices() > simp_percentage * initNumV) )
	int cnt = 0;
	double final = simp_percentage* (double)initNumV;

	std::cout << "heap iteration" << std::endl;
	while ( heap.size() && tot >= final )
	{
		auto tmp = heap.top(); heap.pop();
		// std::cout << tot   <<"|"<< final << " |cost:" << tmp.cost << std::endl;
		auto vhat = tmp.vhat;

		contract_v2(tmp.vi, tmp.vj, vhat);
		cnt++ ;
	}
}



void Mesh::add_pair(int i, int j)
{
	// std::cerr << "add " << i << "|" << j << "\n";
	assert(i != j);
	if (i > j) std::swap(i, j);

	double _d =  vertices[i].cal_dist(vertices[j]);
	// std::cerr << i << "|dist:" << _c << "|" << thre_dist << "\n";
	// if ( _d < thre_dist){

		// cal cost of possible contraction
		// we choose v1+v2 as final Q
		// v hat is the (v1+v2) / 2
		// Eigen::Matrix4d Q = findQ[v1] + findQ[v2];
	Eigen::Vector4d newx;
	// update Q, v, cost
	double cost = cal_cost(i,j, newx);
	// std::cerr << "cost:" << cost << "\n";
	auto vp = ValidPair(i, j, cost, newx );
	heap.push(vp);
	// }
}




Mesh *mesh;
int main(int argc,char *argv[])
{
	// if(argc != 5) std::cout << "must send 4 params" << "\n";

	std::string mesh_path = argv[1];
	std::string output = argv[2];
	simp_percentage = atof(argv[3]);
	t = atof(argv[4]);
	std::cerr << simp_percentage << "\n";

	auto start = std::chrono::high_resolution_clock::now();
	mesh = new Mesh();
	mesh->loadMesh(mesh_path, mesh);   
	
	initNumV = mesh->numVertices();
	toDel.resize(initNumV); std::fill(toDel.begin(), toDel.end(), 0);
	int init_face_cnt = mesh->numFaces();
	//  iters;
	// thre_dist= 1.0;
	// simp_percentage = 0.5, 
	mesh->simp();
    mesh->writeMesh(output, mesh);

	auto end = std::chrono::high_resolution_clock::now();
	std :: cerr << "time comsumption: " <<  (end- start).count() / 1000000 / 10000.0 << std::endl;
	
	std::cerr << " init Vertex cnt: " << initNumV << "\n";
	std::cerr << " init Face cnt: " << init_face_cnt << "\n";
	std::cerr << " post Vertex cnt: " << post_v_cnt << "\n";
	std::cerr << " post Face cnt: " << mesh->numFaces() << "\n";

    return 0;

}