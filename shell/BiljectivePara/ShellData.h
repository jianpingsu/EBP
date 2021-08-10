#ifndef SHELLDATA_H
#define SHELLDATA_H

#include <Eigen/Dense>
#include <map>
#include"Common.h"
#include"TriangleInterface.h"
struct ShellData
{
public:
	ShellData()
	{
		dim = 2;
		mv_num = 0;
		mf_num = 0;
		sf_num = 0;
		sv_num = 0;
		mesh_measure = 0;
		w_uv.resize(0, 0);
	};

	void add_new_patch(const Eigen::MatrixXd&, const Eigen::MatrixXi&, const Eigen::RowVectorXd &center);

	void mesh_improve();

	void update_shell();

	double shell_factor = 10;

	double energy;		// mesh energy

	long mv_num, mf_num;
	long sv_num, sf_num;
	Eigen::MatrixXd m_V;	// input initial mesh V
	Eigen::MatrixXi m_T;	// input initial mesh F/T

	Eigen::MatrixXd w_uv;	// whole domain uv: mesh + free vertices
	Eigen::MatrixXi s_T;	// shell domain tets: shell tets
	Eigen::MatrixXi w_T;

	Eigen::VectorXd m_M;	// mesh area or volume
	Eigen::VectorXd s_M;	// shell area or volume
	Eigen::VectorXd w_M;	// area/volume weights for whole
	Eigen::MatrixXi surface_F;
	double mesh_measure;	// area or volume
	long v_num;
	long f_num;

	Eigen::VectorXi frame_ids;
	Eigen::MatrixXd frame_V;
	Eigen::VectorXi internal_bnd;

	std::vector<int> component_sizes;	// multi-chart support
	std::vector<int> bnd_sizes;

	int dim = 2;						// dimension for ambient space. Same for mesh/shell
};


#endif //SHELLDATA_H
