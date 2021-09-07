#include "PolyMesh/IOManager.h"
#include "Eigen/Dense"
#include <queue>
#include <map>
#include <time.h>

#define lambda 1000

using namespace acamcad;
using namespace polymesh;
using namespace std;

struct Edge_priority{
    double cost;
    MEdge* e;
    MPoint3 new_point;

};
struct cmp {
    bool operator()(const Edge_priority& a, const Edge_priority& b) const { return a.cost < b.cost; }
};


priority_queue<Edge_priority, vector<Edge_priority>, cmp> Cost;
map<MVert*, Eigen::Matrix4d> Q_v;

void cal_Q(MVert* v, PolyMesh* mesh){
    // if (mesh->isBoundary())
    Eigen::Matrix4d Q_ans;
    Q_ans.setZero();

    for ( VertexFaceIter it = mesh->vf_iter(v); it.isValid(); ++it){
        MVector3 normal = (*it)->normal();
        double d = - dot(v->position(), normal);
        Eigen::Matrix<double, 4, 1> p (normal[0], normal[1], normal[2], d);
        Q_ans += p * p.transpose();
    }
    Q_v[v] = Q_ans;
}

void cal_cost(MEdge* e){
    MHalfedge* he = e->halfEdge();
    MVert* from = he->fromVertex(); MVert* to = he->toVertex();
    Eigen::Matrix4d Q_plus = Q_v[from] + Q_v[to];
    
    MPoint3 new_point = from->position() + to->position();
    Eigen::Vector4d new_vec = {new_point[0], new_point[1], new_point[2], 1};
    Edge_priority tmp;
    tmp.cost = new_vec.transpose() * Q_plus * new_vec;
    tmp.new_point = new_point;
    tmp.e = e;
    Cost.push(tmp);
}

void collapse(Edge_priority& ep, PolyMesh* mesh){
    MHalfedge* he = ep.e->halfEdge();
    MHalfedge* he_oppo = he->pair();
	MVert *v_from = he->fromVertex(), *v_to = he->toVertex();
    if (mesh->is_collapse_ok(he)){
		v_to->setPosition(ep.new_point);
        mesh->collapse(he);
    }
    else if (mesh->is_collapse_ok(he_oppo)){
		v_from->setPosition(ep.new_point);
        mesh->collapse(he_oppo);
    }
}

int main(int argc, const char **argv){

    if (argc != 3){
        cout << "./a i.obj o.obj" << endl;
        return -1;
    }
    PolyMesh* mesh = new PolyMesh();
    string input = argv[1];
    string output = argv[2];
    loadMesh(input, mesh);
    int init_num_faces = mesh->numPolygons();
    int ratio = 0.5;

    // calculate Q for each vertice
    for (MVert* v: mesh->vertices()){
        cal_Q(v, mesh);
    }
    for (MEdge* e: mesh->edges()){
        cal_cost(e);
    }

    while (mesh->numPolygons() > init_num_faces * ratio){
        // select valid pairs
        Edge_priority e = Cost.top(); Cost.pop();
        cout << e.cost << endl;
        collapse(e, mesh);
    }

    writeMesh(output, mesh);
    return 0;
}


