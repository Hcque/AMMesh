#pragma once
#include <PolyMesh/PolyMesh.h>
#include <PolyMesh/PolyMesh_Base.h>
#include <vector>
using namespace acamcad;
using namespace polymesh;

struct pathInfo{
	int from, to; // in completde g
	double dist; // weight of compledted g
	int selected;
	std::vector<int> path; 
	bool operator<(const pathInfo& b){
		return dist < b.dist;
	}
	pathInfo(){}
	pathInfo(int _f, int _t, double _d, std::vector<int>& _p): 
		from(_f), to(_t), dist(_d), path(_p), selected(0) {}
};

//��������lmk֮��ĵ��·��
void Dijkstra_group(PolyMesh& Mesh, std::vector<int>& lmk, std::vector<std::vector<int>>&path);

//��������lmk֮��ĵ��·��
void Dijkstra(PolyMesh& Mesh, int s_p, int e_p, std::vector<int>& path);