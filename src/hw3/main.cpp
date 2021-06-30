
#include <iostream>
#include <fstream>
#include <sstream>
#include "PolyMesh/IOManager.h"
#include <string>

using namespace acamcad;
using namespace polymesh;
using namespace std;

int main(int argc, char** argv)
{

    if (argc != 4){
        cout << "./exe .obj 1 8" << endl;
        return -1;
    }
    string mesh_path = argv[1];
    int sigma = atoi(argv[2]);
    int iter_num = atoi(argv[3]);

    PolyMesh *mesh = new PolyMesh();
    loadMesh(mesh_path, mesh);
    mesh->updateMeshNormal();

    vector<MVector3> new_normal(mesh->numPolygons());
    vector<double> face_area(mesh->numPolygons());
    vector<MPoint3> face_center(mesh->numPolygons());

    for (MPolyFace *f: mesh->polyfaces() ){
        // cout << f << endl;
        int index = (*f).index();
        new_normal[index] = f->normal();
        face_center[index] = mesh->calculatFaceCenter(f);

        vector<MVert*> points3;
        for (FaceVertexIter fvit = mesh->fv_iter(f); fvit.isValid(); ++fvit){
            points3.push_back((*fvit));
        }
        auto e20 =  points3[2]->position() - points3[0]->position();
        auto e10 =  points3[1]->position() - points3[0]->position();
        face_area[index] = cross(e20,e10).norm() * 0.5 ;
    }

    double sigma_c = 0.0;
    for (MPolyFace *f: mesh->polyfaces() ){
        // cout << f << endl;
        int index = (*f).index();
        for (FaceFaceIter ff = mesh->ff_iter(f); ff.isValid(); ++ff){
            sigma_c += (face_center[index] - face_center[(*ff)->index()]).norm();
            // cout << "sigma_c:" << sigma_c << endl;
        }
    }
    sigma_c /= mesh->numPolygons() * 3 * 2;
    cout << "sigma_c:" << sigma_c << endl;

    // update normal
    for (MPolyFace *f: mesh->polyfaces() ){
        int index_i = f->index();
        double Ws = 0;
        double Wr = 0;
        double Aj = 0;
        double Kp = 0;
        MVector3 ans(0.0, 0.0, 0.0);
        sigma_c  = 1;

        for (FaceFaceIter ff = mesh->ff_iter(f); ff.isValid(); ++ff){
            int index_j = (*ff)->index();
            double d_center = (face_center[index_i] - face_center[index_j]).norm();
            double d_norm = (new_normal[index_i] - new_normal[index_j]).norm();
            // cout << face_center[index_i] << face_center[index_j] << endl;
            // cout << d_center << d_norm << endl;
            Ws = exp((-1.0)* (d_center * d_center) / (2*sigma_c*sigma_c));
            Wr = exp((-1.0)* (d_norm * d_norm) / (2*sigma*sigma));
            Aj = face_area[index_j];
            Kp += Wr*Ws*Aj;
            // cout <<  Ws << Wr << Aj << Kp << endl;
            ans += Aj*Wr*Ws * new_normal[index_j];
        }
        new_normal[index_i] = (ans / Kp).normalized();
        // cout << ans[0] << endl;
    }

    // mesh->vertices
    // update vertex
    for (int i = 0; i < iter_num; i++){
        for (MVert* v: mesh->vertices()){
            int index_i = v->index();
            MPoint3 xi = v->position();
            cout << v->x() << " " ;

            MPoint3 delta(0.0, 0.0, 0.0);
            int n = 0;
            for (VertexFaceIter vfit = mesh->vf_iter(v); vfit.isValid(); ++vfit){
                n++;
                int index_j = (*vfit)->index();
                MVector3 nj = new_normal[index_j];
                MPoint3 cj = mesh->calculatFaceCenter(*vfit);
                delta = delta + dot(nj, (cj - xi)) * nj.point();
            }
            v->setPosition(xi + delta/(2*n*n));
            cout << v->x() << endl;
        }
    }

    writeMesh("outpu.obj", mesh);
    return 0;
    
}