
#include <iostream>
#include <fstream>
#include <sstream>
#include "PolyMesh/IOManager.h"
#include <string>
#include <vector>

using namespace acamcad;
using namespace polymesh;
using namespace std;
int T = 20;

std::vector<vector<MVert*>> vers;
std::vector<MVector3> normals;

inline double dist(const MPoint3& v1, const MPoint3& v2){
    return v1.distance(v2);
}
void test_(PolyMesh *mesh){
    const int F = mesh->numPolygons();
    vers.resize(F);

    for (int i=0;i<mesh->numPolygons();i++){
        MPolyFace* f = mesh->polyfaces()[i];
        int faceid = f->index();

        for (auto it = mesh->fv_iter(f); it.isValid();++it){
            (*it)->position();
        }
        auto e = f->halfEdge();
        vers[faceid].push_back(e->fromVertex());
        e=e->next(); vers[faceid].push_back( e->fromVertex());
        e=e->next(); vers[faceid].push_back( e->fromVertex());

     
    }

    while (T--){
        std::cerr<<"T"<<T<<"\n";

    for (int i=0;i< mesh->numVertices();i++){
        auto v = mesh->vertices()[i];

        MVector3 N(0,0,0); double c = 0.0;
        for (auto vf = mesh->vf_iter(v);vf.isValid();++vf){
            N += (*vf)->normal(); c+=1;
        }
        N /= c;
        N.normalize();

        MVector3 dv(0,0,0); c = 0;
        for (auto vf = mesh->vf_iter(v);vf.isValid();++vf){
            for (auto fv = mesh->fv_iter((*vf));fv.isValid();++fv){
                MVector3 h = (*fv)->position() - v->position();
            // if (h.dot(h) < 0.5){
                dv += h.dot(N); c+=1;
            }
        }


        v->setPosition(v->position()+dv/c);
    }

    }

	std::cerr << "Test done\n" ;
}





int amain(int argc, char** argv)
{
    if (argc < 2)
	{
		std::cout << "========== Hw3 Usage  ==========\n";
		std::cout << std::endl;
		std::cout << "Input:	ACAM_mesh_HW2.exe	mesh.obj\n";
		std::cout << std::endl;
		std::cout << "=================================================\n";
		return -1;
	}

    std::string mesh_path = argv[1];
	PolyMesh* mesh = new PolyMesh();
	loadMesh(mesh_path, mesh);
	std::cout << "num verts: " << mesh->numVertices() << std::endl;

	test_(mesh);

    writeMesh("./result.obj", mesh);
    // test2_(mesh);
	delete mesh;
    return 0;
    
}