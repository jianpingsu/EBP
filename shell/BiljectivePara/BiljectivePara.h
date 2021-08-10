#pragma once
#include"Parafun.h"
#include"ShellData.h"

#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshDefinition.h"

using namespace Eigen;
using namespace std;

class BiljectivePara
{
public:
	
	BiljectivePara(string filename);
	~BiljectivePara();

	void parameterization();
	void load();
	void write_obj();
	void shelltri(MatrixXi& tri, MatrixXd& pos);
	double adjust_weight(double conv_mesh, double last_energy);

protected:

	string		modelname;
	Mesh		mesh;
	ShellData	shell_data;
	std::shared_ptr<Parafun> parafun_solver = nullptr;

	double convgence_con_rate = 1e-5;
	int MAX_ITER_NUM = 5000;
	bool is_initialization;

	bool weight1;
	bool weight2;
};