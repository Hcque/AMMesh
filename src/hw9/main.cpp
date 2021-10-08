
// caution summary:

// set <opertoar
// insert to set , const &
// local for faces
// idx needs repad when writing
// costomized heap, Hash


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

using namespace acamcad;
using namespace polymesh;


int iters;
double thre_dist;
double simp_percentage ;
int initNumV;


struct ValidPair
{
	// MEdge* me;
	int vi, vj;
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

// https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key

  struct HashFunction
  {
    std::size_t operator()(const ValidPair& k) const
    {
      using std::size_t;
      using std::hash;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:

      return (hash<int>()(k.vi) ^ (hash<int>()(k.vj) << 1) );
    }
  };

};


std::unordered_map<ValidPair, int, ValidPair::HashFunction> pairDup;
std::vector<int> toDel;

// template <class T>
class Heap : public std::priority_queue<ValidPair>
{
public:
	bool remove_related(const int v)
	{
		bool flag = false;
		for ( auto it = this->c.begin(); it != this->c.end(); ++it)
		{
			if ( (*it).vi == v || (*it).vj == v ){
				this->c.erase(it);
				flag = true;
			}	
		}
	
		std::make_heap(this->c.begin(), this->c.end(), this->comp);
		return flag;
	}
};

std::unordered_map<int , Eigen::Matrix4d> findQ;
Heap heap;


class Mesh
{
public:

	struct Face{
		mutable int dim[3];
		Face(int _i, int _j, int _k){
			dim[0] = _i;
			dim[1] = _j;
			dim[2] = _k;
			sort();
		}
		void sort(){
			std::sort(dim,dim+3);
		}

		MVector3 normal() const
		{
			MVector3 ba = dim[1]-dim[0];
			MVector3 ca = dim[2]-dim[0];
			ba.normalize(); ca.normalize();
			return ba.cross(ca);
		}
		bool count(int vi) const { 
			return std::find(dim,dim+3, vi) != dim+3;
		}

		int find_other(int v1, int v2) const
		{
			std::set<int> copy(std::begin(dim), std::end(dim));
			copy.erase(v1);
			copy.erase(v2);
			return *copy.begin();
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
			if (dim[0] == v1) o2 == dim[1], o3 = dim[2];
			else if (dim[1] == v1) o2 = dim[0], o3 = dim[2];
			else if (dim[2] == v1) o2 = dim[0], o3 = dim[1];
			// else assert(0);
			else return false;
			return true;
		}
		auto getEdges()const 
		{
			std::vector<ValidPair> ans;
			ans.push_back(ValidPair(dim[0], dim[1]));
			ans.push_back(ValidPair(dim[1], dim[2]));
			ans.push_back(ValidPair(dim[0], dim[2]));
			return ans;
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
		int index;
		std:: set<Face> fas;
		Vertice(){}
		Vertice(int i,int j, int k, int idx){
			position[0] = i;
			position[1] = j;
			position[2] = k;
			index = idx;
		}

		void add_face(const Face& f){ fas.insert(f); }
		void remove_face(const Face& f){ fas.erase(f); }

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

public:
	void loadMesh(std::string, Mesh* mesh);
	void writeMesh(std::string, Mesh* mesh);
	int numVertices() const { return vertices.size(); }
	void contract_v2(int v1, int v2, Eigen::Vector4d);
	
	std::vector<Face> faces;
	std::vector<Vertice> vertices;

	void addFace(int i, int a,int b,int c)
	{
		vertices[i].fas.insert(Face(a,b,c));
	}

	void addVertice(int i, int j, int k, int index)
	{
		vertices.push_back(Vertice(i,j,k,index));
	}

	// void update_Node(int vid, MPoint3 newX)
	// {
	// 	vertices[vid].position = newX;
	// }
};


void Mesh::contract_v2(int v1, int v2, Eigen::Vector4d vhat)
	{
		auto faceset = vertices[v2].fas; // f1,f2, 5,6,7
		for (auto &ff : faceset) {
			if (ff.count(v1)) { // f1,f2
				int vother = ff.find_other(v1,v2);
				vertices[vother].remove_face(ff);
				vertices[v2].remove_face(ff);
				vertices[v1].remove_face(ff);
			}
			else{  // 5,6,7
			assert(ff.count(v2));
				if (!ff.count(v2)){
					// vertices[v2].fas
					std::cout << v2 << "\n";
					for (auto & f: vertices[v2].fas) std::cout << f <<" ";
					std::cout << std::endl;
				}

				int o2,o3;
				ff.find_other2(v2,o2,o3);

				Face newf(v1,o2,o3);
				vertices[v1].add_face(newf);
				
			}
		}
		std::cerr << "constraction" << v2 << "\n";
		// vertices
		toDel[v2] = 1;
		// connect neigbours  to v1
		for (int i: vertices[v1].oneRingNeibours() )
			vertices[i].reset(v2,v1);
		// update v1
		vertices[v1].position[0] = vhat[0];
		vertices[v1].position[1] = vhat[1];
		vertices[v1].position[2] = vhat[2];
	}


int init_face_cnt = 0;
Mesh *mesh;
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
				mesh->addVertice(a,b,c, ++cnt);
				break;
			case 'f':
				ss >> fa >> fb >> fc;
				mesh->addFace(fa-1, fa-1,fb-1,fc-1);
				mesh->addFace(fb-1, fb-1,fc-1,fa-1);
				mesh->addFace(fc-1, fc-1,fa-1,fb-1);
				init_face_cnt ++;
				break;
		}
	}
}

std::vector<int> mapidx;
std::unordered_set<Mesh::Face, Mesh::Face::HashFunction> fascnt;
int post_v_cnt = 0;
int post_f_cnt = 0;
void Mesh::writeMesh(std::string mesh_path, Mesh* mesh)
{
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
			fascnt.insert(ff);

					out << 'f' << " " <<
			mapidx[ff.dim[0]]+1 << " " <<
			mapidx[ff.dim[1]]+1 << " " <<
			mapidx[ff.dim[2]]+1 << std::endl;
			post_f_cnt ++ ; 
		}
	}
}



Eigen::Matrix4d cal_Q(int vid)
{
	Eigen::Matrix4d ans;

	MPoint3 X = mesh->vertices[vid].position;

	for (auto &ff : mesh->vertices[vid].fas)
	{
		MVector3 normal  = ff.normal();
		assert( normal.dot(normal) - 1.0 < 0.001);

		double d = -1.0f * ( normal[0]*X[0] + normal[1]*X[1] + normal[2]*X[2] );
		Eigen::Matrix<double, 4,1> p(normal[0], normal[1], normal[2], d);
		Eigen::Matrix4d Qi = p * p.transpose();
		ans += Qi;
	}
	return ans;
}

double cal_dist(int v1, int v2)
{

	MPoint3 X = mesh->vertices[v1].position;
	MPoint3 Y = mesh->vertices[v2].position;

	return X.distance(Y);
	return 0.0;
}

double cal_cost(int v1, int v2, const Eigen::Matrix4d& Q, Eigen::Vector4d& newx )
{
	MPoint3 X = mesh->vertices[v1].position;
	MPoint3 Y = mesh->vertices[v2].position;

	Eigen::Vector4d X4(X[0],X[1], X[2], 1.0f);
	Eigen::Vector4d Y4(Y[0],Y[1], Y[2], 1.0f);
	Eigen::Vector4d xhat = (X4+Y4) / 2;
	newx = xhat;
	return std::fabs( xhat.transpose() * Q * xhat );
}



void simp(Mesh* mesh)
{
	// cal Q, stored in map
	std::cout << "cal Q, stored in map" << std::endl;

	for (int i=0;i< mesh->vertices.size(); i++) {
		int vid = i;

		Eigen::Matrix4d Q = cal_Q(vid);
		findQ[vid] = Q;
		// std::cout << vid << "|" << Q.data << "\n";
	}
	// collect valid pairs, push in ti heap;
	std::cout << "ccollect valid pairs, push in ti heapp" << std::endl;

	for (auto v: mesh->vertices) for (auto &f : v.fas){
		for (auto& edge: f.getEdges()){

			if (pairDup[edge] == 1) continue;
			pairDup[edge] = 1;
			int v1 = edge.vi, v2 = edge.vj;

			double _c =  cal_dist(v1,v2);
			// std::cerr << "cost:" << _c << "\n";
			if ( _c < thre_dist){

				// cal cost of possible contraction
				// we choose v1+v2 as final Q
				// v hat is the (v1+v2) / 2
				Eigen::Matrix4d Q = findQ[v1] + findQ[v2];
				Eigen::Vector4d newx;
				double cost = cal_cost(v1,v2,Q, newx);
				std::cerr << "cost:" << cost << "\n";

				heap.push(ValidPair(v1, v2, cost, newx ));
				
			}
		// std::cout << heap.size() << " ";
		}
	}
		std::cout <<  "heap size: " << heap.size() << "\n";

	// iters
	// while ( heap.size() && ( mesh->numVertices() > simp_percentage * initNumV) )
	int cnt = 0;

	std::cout << "cheap iter" << std::endl;
	while ( heap.size() && (cnt++ < iters) )
	{
		auto tmp = heap.top(); heap.pop();

		// std::cout << cnt << " |cost:" << tmp.cost << std::endl;
		auto vhat = tmp.vhat;

		mesh->contract_v2(tmp.vi, tmp.vj, vhat);
		
		// remove pairs linked to v1 or v2;
		heap.remove_related(tmp.vi);		
		heap.remove_related(tmp.vj);		


		// to be continued
		// cal Q of v hat
		// updateQ();

		// add possible valid pairs corrosponding to hhat


	}
}



int main(int argc,char *argv[])
{
	// if(argc != 5) std::cout << "must send 4 params" << "\n";

	std::string mesh_path = argv[1];
	std::string output = argv[2];
	iters = atoi(argv[3]);
	thre_dist = atoi(argv[4]);

	mesh = new Mesh();
	mesh->loadMesh(mesh_path, mesh);   
	
	int n = mesh->numVertices();
	std::cout << "num verts: " << n << std::endl;
	toDel.resize(n); std::fill(toDel.begin(), toDel.end(), 0);

	//  iters;
	// thre_dist= 1.0;
	// simp_percentage = 0.5 , 
	initNumV = n;

	simp(mesh);

    mesh->writeMesh("./result.obj", mesh);

    // test2_(mesh);
	delete mesh;


	std::cerr << " init Vertex cnt: " << initNumV << "\n";
	std::cerr << " init Face cnt: " << init_face_cnt << "\n";
	std::cerr << " post Vertex cnt: " << post_v_cnt << "\n";
	std::cerr << " post Face cnt: " << post_f_cnt << "\n";

    return 0;

}