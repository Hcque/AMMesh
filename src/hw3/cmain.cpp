
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

std::vector<MVert*> v1,v2,v3;
std::vector<MVert*> vers[3];
std::vector<MVector3> normals;

inline double dist(const MPoint3& v1, const MPoint3& v2){
    return v1.distance(v2);
}
void test_(PolyMesh *mesh){

    for (int i=0;i<mesh->numPolygons();i++){
        auto f = mesh->polyfaces()[i];
        auto e = f->halfEdge();
        v1.push_back(e->fromVertex());
        e=e->next(); v2.push_back( e->fromVertex());
        e=e->next(); v3.push_back( e->fromVertex());

        // std::cerr<<v2.size() << "\n";
        MVector3 ans(0,0,0); double c=0;
        for (auto it = mesh->ff_iter(f); it.isValid();it++, c+=1){
            ans += it.cur_pointer()->normal();
        }
        ans /= c;
		// std::cerr << ans[0] << "|" << ans[1] << "|" << ans[2] << "\n";
        normals.push_back(ans);
    }
    vers[0] = (v1);
    vers[1] = (v2);
    vers[2] = (v3);

    while (T--){
        std::cerr<<"T"<<T<<"\n";

    for (int i=0;i< mesh->numVertices();i++){
        auto v = mesh->vertices()[i];
        MVector3 dv(0,0,0); double c = 0;
        // closet normal
        MVector3 N; double distance=1e9;
        for (auto vf = mesh->vf_iter(v);vf.isValid();++vf){
            double cur = dist((*vf)->getFaceCenter(),v->position());
            if (cur < distance) {
                distance = cur;
                N = (*vf)->normal();
            }
        }

        for (auto vv = mesh->vv_iter(v);vv.isValid();++vv){
            MVector3 h = (*vv)->position() - v->position();
            // if (h.dot(h) < 0.5){
                dv += h.dot(N); c+=1;
        }

        v->setPosition(v->position()+dv/c);
    }

    }

	std::cerr << "Test done\n" ;
}





int cmain(int argc, char** argv)
{
    /*
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
    test2_(mesh);
	delete mesh;
    */

   	if (argc != 5)
	{
		std::cout << "========== Hw3 Usage  ==========\n";
		std::cout << std::endl;
		std::cout << "ACAM_mesh_HW3.exe [model] [SigmaNormal] [NormalIterNum] [VertexIterNum]\n";
		std::cout << "ACAM_mesh_HW3.exe	mesh.obj 1 8 8\n";
		std::cout << std::endl;
		std::cout << "=================================================\n";
		return -1;
	}
	//������
	std::string mesh_path = argv[1];
	std::stringstream ss;
	std::string cmd2 = argv[2], cmd3 = argv[3], cmd4 = argv[4];
	ss << cmd2 + " " + cmd3 + " " + cmd4;
	double SigmaCenter, SigmaNormal;
	double NormalIterNum, VertexIterNum;
	ss >> SigmaNormal >> NormalIterNum >> VertexIterNum;

	PolyMesh* mesh = new PolyMesh();
	loadMesh(mesh_path, mesh);

	// mesh->updateMeshNormal();
	std::vector<MVector3> NewNormal(mesh->numPolygons());//ÿ����ķ���
	std::vector<double> FaceArea(mesh->numPolygons());//ÿ��������
	std::vector<MPoint3> FaceCenter(mesh->numPolygons());//

	for (MPolyFace* fh : mesh->polyfaces())
	{
		int f_id = (*fh).index();
		NewNormal[f_id] = (*fh).normal();
		std::vector<MVert*> P;
		for (FaceVertexIter vv_it = mesh->fv_iter(fh); vv_it.isValid(); ++vv_it)
		{
			P.push_back(*vv_it);
		}
		auto e12 = P[1]->position() - P[0]->position();
		auto e13 = P[2]->position() - P[0]->position();
		double area = cross(e12, e13).norm() * 0.5;
		FaceArea[f_id] = area;
		FaceCenter[f_id] = mesh->calculatFaceCenter(fh);
	}

	SigmaCenter = 0;
	for (MPolyFace* fh : mesh->polyfaces())
	{
		int f_id = (*fh).index();
		for (FaceFaceIter nei_fh = mesh->ff_iter(fh); nei_fh.isValid(); ++nei_fh)
		{
			int ff_id = (*nei_fh)->index();
			SigmaCenter += (FaceCenter[f_id] - FaceCenter[ff_id]).norm();
		}
	}
	SigmaCenter /= mesh->numPolygons() * 3;

	for (int i = 0; i < NormalIterNum; i++)
	{
		for (MPolyFace* fh : mesh->polyfaces())
		{
			double Kp = 0;
			MVector3 NewN(0, 0, 0);
			int fh_id = (*fh).index();
			for (FaceFaceIter nei_fh = mesh->ff_iter(fh); nei_fh.isValid(); ++nei_fh)
			{
				int nei_fh_id = (*nei_fh)->index();
				double delta_center = (FaceCenter[fh_id] - FaceCenter[nei_fh_id]).norm();
				double delta_normal = (NewNormal[fh_id] - NewNormal[nei_fh_id]).norm();
				double Aj = FaceArea[nei_fh_id];
				double Ws = exp(-delta_center * delta_center / (2 * SigmaCenter * SigmaCenter));
				double Wr = exp(-delta_normal * delta_normal / (2 * SigmaNormal * SigmaNormal));
				NewN += Aj * Ws * Wr * NewNormal[nei_fh_id];
				Kp += Aj * Ws * Wr;
			}
			NewNormal[fh_id] = NewN / Kp;
			NewNormal[fh_id] /= NewNormal[fh_id].norm();
		}
	}

	for (int i = 0; i < VertexIterNum; i++)
	{
		for (MVert* vh : mesh->vertices())
		{
			MPoint3 x_i = (*vh).position();
			MPoint3 delta_xi(0, 0, 0);
			int Nei_count = 0;
			for (VertexFaceIter fh = mesh->vf_iter(vh); fh.isValid(); ++fh)
			{
				Nei_count++;
				MPoint3 cj = mesh->calculatFaceCenter(*fh);
				MVector3 nj = NewNormal[(*fh)->index()];
				delta_xi = delta_xi + nj * (nj.data()[0] * (cj - x_i).data()[0] + nj.data()[1] * (cj - x_i).data()[1] + nj.data()[2] * (cj - x_i).data()[2]);
			}
			x_i = x_i + delta_xi / Nei_count;
			(*vh).setPosition(x_i);
		}
	}

	cout << "output result" << endl;
	writeMesh("result.obj", mesh);
	return 0;

   
    
}