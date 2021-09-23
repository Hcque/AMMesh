#include <iostream>
#include <fstream>
#include <sstream>
#include "PolyMesh/IOManager.h"
#include <string>
#include <cmath>

using namespace acamcad;
using namespace polymesh;
using namespace std;

double tuo;
inline double W(double x)
{
	// double tuo = 0.1;
	return exp(-(x*x)/(2* tuo*tuo));
}

int main(int argc, char** argv)
{
	if (argc != 5)
	{
		std::cout << "========== Hw3 Usage  ==========\n";
		std::cout << std::endl;
		std::cout << "ACAM_mesh_HW3.exe [model] [SigmaNormal] [NormalIterNum] [VertexIterNum]\n";
		std::cout << "ACAM_mesh_HW3.exe	mesh.obj iter res.obj\n";
		std::cout << std::endl;
		std::cout << "=================================================\n";
		return -1;
	}
	//������
	std::string mesh_path = argv[1];
	std::stringstream ss;
	std::string cmd2 = argv[2];
	std::string cmd3 = argv[3];
	std::string cmd4 = argv[4];
	double SigmaCenter, SigmaNormal;
	double NormalIterNum, VertexIterNum;
	int iter = atoi( cmd2.c_str());
	tuo = atof(cmd4.c_str());

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


	while (iter -- ){
			// cout << "iter:" << iter << endl;
		for (MVert* vh : mesh->vertices())
		{
			MPoint3 x_i = (*vh).position();

			// avg normal for that vertex
			int Nei_count = 0;
			MVector3 avg_normal(0,0,0);
			for (VertexFaceIter fh = mesh->vf_iter(vh); fh.isValid(); ++fh)
			{
				Nei_count++;
				avg_normal += (*fh)->normal();
			}
			avg_normal /= Nei_count;
			avg_normal.normalize();

			// d
			MVector3 sum_(0,0,0);
			double normalize = 0;
			for (VertexVertexIter vv = mesh->vv_iter(vh); vv.isValid(); ++ vv)
			{
				MVector3 dpoint = (*vv)->position() - x_i;
				double hh = dpoint.dot(avg_normal);
				double w1 = W( sqrt(dpoint.dot(dpoint)) );
				double w2 = W( hh );
				sum_ += w1*w2*hh*avg_normal;
				normalize += 1;	
			}
			MVector3 var = sum_ / normalize;
			cout << var[0] << endl;
			(*vh).setPosition( (*vh).position() + var );
		}
	}

	cout << "output result" << endl;
	writeMesh(cmd3, mesh);
	return 0;
}