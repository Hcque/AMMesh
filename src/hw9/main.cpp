#include "PolyMesh/IOManager.h"
#include "Eigen/Dense"
#include <queue>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <time.h>

using namespace acamcad;
using namespace polymesh;


int iters;
double thre_dist;
double simp_percentage ;
int initNumV;

class Vertice
{
public:
	MPoint3 position;
	int index;
	Vertice() {}
	Vertice(int _a, int _b, int _c){
		position[0] = _a;
		position[1] = _b;
		position[2] = _c;
	}
};

struct Face{
	int i,j,k;

};
class Mesh
{
public:
	void loadMesh(std::string, Mesh* mesh);
	void write_mesh(std::string, Mesh* mesh);
	int numVertices() const { return vertices.size(); }

	std::vector<Vertice> vertices;
	std::vector<Face> faces;

};

Mesh *mesh;
void Mesh::loadMesh(std::string mesh_path, Mesh* mesh)
{
	double a, b, c; int fa,fb,fc;
	char ch;

	std::cout << "load..." << "\n";
	std::ifstream in(mesh_path);

	while (in >> ch){

		// printf("%c\n", ch);
		switch (ch)
		{
			case 'v':
				in >> a >> b >> c;
				// std::cout << a << b << c << "\n";
				vertices.push_back(Vertice(a,b,c));
				break;
			case 'f':
				in >> fa >> fb >> fc;
				faces.push_back({fa,fb,fc});
				break;
		}

	}
}

struct ValidPair
{
	// MEdge* me;
	int vi, vj;
	double cost;
	Eigen::Vector4d vhat;

	bool operator<(const ValidPair& other) const { return cost > other.cost; } 
};


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
/*

Eigen::Matrix4d cal_Q(int vid)
{
	Eigen::Matrix4d ans;

	MPoint3 X = mesh->vertices()[vid]->position();
    // neibour faces
	for (auto it = mesh->vf_iter(mesh->vertices()[vid]); it.isValid(); ++it)
	{
		MVector3 normal  = (*it)->normal();
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

	MPoint3 X = mesh->vertices()[v1]->position();
	MPoint3 Y = mesh->vertices()[v2]->position();

	return X.distance(Y);
	return 0.0;
}

double cal_cost(int v1, int v2, const Eigen::Matrix4d& Q, Eigen::Vector4d& newx )
{
	MPoint3 X = mesh->vertices()[v1]->position();
	MPoint3 Y = mesh->vertices()[v2]->position();

	Eigen::Vector4d X4(X[0],X[1], X[2], 1.0f);
	Eigen::Vector4d Y4(Y[0],Y[1], Y[2], 1.0f);
	Eigen::Vector4d xhat = (X4+Y4) / 2;
	newx = xhat;
	return xhat.transpose() * Q * xhat;
}

void simp(PolyMesh* mesh)
{

	// cal Q, stored in map

	std::cout << "cal Q, stored in map" << std::endl;

	for (MVert* mv : mesh->vertices()){
		int vid = mv->index();

		Eigen::Matrix4d Q = cal_Q(vid);
		findQ[vid] = Q;
	}
	// collect valid pairs, push in ti heap;
	std::cout << "ccollect valid pairs, push in ti heapp" << std::endl;

	for (MEdge* me : mesh->edges()){

		if (mesh->isBoundary(me->halfEdge())) 
		{
			std::cout << "boundary" << "\n";
			continue;
		}
		auto v1 =  me->getVert(0)->index();
		auto v2 =  me->getVert(1)->index();

		if ( cal_dist(v1,v2) < thre_dist){

			// cal cost of possible contraction
			// we choose v1+v2 as final Q
			// v hat is the (v1+v2) / 2
			Eigen::Matrix4d Q = findQ[v1] + findQ[v2];
			Eigen::Vector4d newx;
			double cost = cal_cost(v1,v2,Q, newx);

			heap.push({me, v1, v2, cost, newx});
		}
		// std::cout << heap.size() << " ";
	}

	// iters
	// while ( heap.size() && ( mesh->numVertices() > simp_percentage * initNumV) )
	int cnt = 0;

	std::cout << "cheap iter" << std::endl;
	while ( heap.size() && (cnt++ < iters) )
	{
		auto tmp = heap.top(); heap.pop();
		auto e = tmp.me;

		std::cout << cnt << " |cost:" << tmp.cost << std::endl;
		auto vhat = tmp.vhat;
		// vhat contract v1,v2;
		MHalfedge* he = e->halfEdge(), *he_neg = he->pair();
		MVert* v1 = he->fromVertex(), *v2 = he->toVertex();


		// remove pairs linked to v1 or v2;
		heap.remove_related(tmp.vi);		
		heap.remove_related(tmp.vj);		


		// to be continued
		// cal Q of v hat
		// updateQ();

		// add possible valid pairs corrosponding to hhat


	}

}
*/



int main(int argc,char *argv[])
{
	if(argc != 5) std::cout << "must send 4 params" << "\n";

	std::string mesh_path = argv[1];
	std::string output = argv[2];
	iters = atoi(argv[3]);
	thre_dist = atoi(argv[4]);

	mesh = new Mesh();
	mesh->loadMesh(mesh_path, mesh);   
	
	int n = mesh->numVertices();
	std::cout << "num verts: " << n << std::endl;

	//  iters;
	thre_dist= 1.0;
	simp_percentage = 0.5 , initNumV = n;

	// simp(mesh);

    // writeMesh(output, mesh);

	std::cout << "after processing, numV: " << mesh->numVertices() << "\n"; 
    // test2_(mesh);
	delete mesh;
    return 0;

}