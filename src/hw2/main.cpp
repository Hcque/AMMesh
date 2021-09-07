#include "PolyMesh/IOManager.h"
#include <iostream>

using namespace acamcad;
using namespace polymesh;

void test_(PolyMesh *mesh){
	for (auto f: mesh->polyfaces()){
		float r = f->normal()[0];
		float g = f->normal()[0];
		float b = f->normal()[0];
		f->setColor(0.4,100,200); 
		f->setVisibility(false);
		std::cerr << r << "|" << g << "|" << b << "\n";

	}
	std::cerr << "Test done\n" ;
}

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		std::cout << "========== Hw2 Usage  ==========\n";
		std::cout << std::endl;
		std::cout << "Input:	ACAM_mesh_HW2.exe	mesh.obj\n";
		std::cout << std::endl;
		std::cout << "=================================================\n";
		return -1;
	}

	//��������
	std::string mesh_path = argv[1];
	PolyMesh* mesh = new PolyMesh();
	loadMesh(mesh_path, mesh);
	std::cout << "num verts: " << mesh->numVertices() << std::endl;

	// mean_Curvature(mesh);
	// abs_mean_Curvature(mesh);
	// gaussian_Curvature(mesh);
	test_(mesh);
	std::string save_path = argv[2];
	writeMesh(save_path, mesh);
	delete mesh;
	return 0;
}