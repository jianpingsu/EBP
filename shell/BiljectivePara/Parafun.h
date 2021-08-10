#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "ShellData.h"
#include "PardisoSolver.h"
#include <time.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshDefinition.h"

using namespace Eigen;
using namespace std;

#define PI 3.141592653589793238463L

class Parafun
{

public:
	Parafun(ShellData & data) :d_(data)
	{
		pardiso = NULL;
		is_first = true;
		bound_distortion_K = 250;
		barrer_coef = d_.mesh_measure * 1e-8;
	};
	~Parafun();

	void after_mesh_improve();
	void init();
	void handle_mintri();
	void init_area();
	void setvirtualtri();
	void Pre_calculate();

	double BPE(bool is_ip_convrate, bool is_slim_convrate);
	void CM(bool is_interp = false);
	void SLIM(bool is_interp = false);
	void Update_source_same_t();

	void fungrid(const VectorXd &x);

	void Energy(const VectorXd &x, double &energy, bool is_interp = false, bool is_whole = true);
	void Energysource();
	double compute_energy(const Eigen::MatrixXd & x, bool whole = false);
	void adjust_shell_weight(double new_weight);

	void max_step(const VectorXd &xx, const VectorXd &dd, double &step);
	void tmaxdetect(const VectorXd &x, const VectorXd &d, double& tmax);
	void backtracking_line_search(const VectorXd &x, const VectorXd &d, const VectorXd &negetive_grad, double &alpha, bool is_interp = false);
	double get_smallest_pos_quad_zero(double a, double b, double c);

	double distance(double s0, double s1, double e0, double e1, double p0, double p1);

	double newton_equation(const double & a, const double & b, const double & K);
	void local_coordinate_inverse(int i, double &p00, double &p01, double &p10, double &p11);
	void local_coordinate_inverse_scaf(int i, double &p00, double &p01, double &p10, double &p11);

	ShellData &d_;
	int	total_num;
	int F_N;
	int V_N;
	int dim = 2;
	vector<int> F0;
	vector<int> F1;
	vector<int> F2;
	VectorXd position_of_mesh;
	vector<double> area;
	vector<double> area_src;
	double area_threshold;

	bool is_first;
	double Intp_T_Min;
	double changetocm_flag;
	double bound_distortion_K;
	vector<double> source_p00;
	vector<double> source_p01;
	vector<double> source_p10;
	vector<double> source_p11;
	vector<double> update_p00;
	vector<double> update_p01;
	vector<double> update_p10;
	vector<double> update_p11;

	PardisoSolver* pardiso;
	vector<int> pardiso_ia;
	vector<int> pardiso_ja;
	vector<double> pardiso_a;
	vector<double> pardiso_b;

	double energy_mesh;
	double energy_barrier;
	double energy_all;
	double energy_shell;

	vector<int> id_h00; vector<int> id_h01; vector<int> id_h02; vector<int> id_h03; vector<int> id_h04; vector<int> id_h05;
	vector<int> id_h11; vector<int> id_h12; vector<int> id_h13; vector<int> id_h14; vector<int> id_h15;
	vector<int> id_h22; vector<int> id_h23; vector<int> id_h24; vector<int> id_h25;
	vector<int> id_h33; vector<int> id_h34; vector<int> id_h35;
	vector<int> id_h44; vector<int> id_h45;
	vector<int> id_h55;

	double barrer_coef;
	double threhold;
	double average_length;

	int BE_N;
	int V_F_N;
	int AV_F_N;
	vector<int>	V_F0;
	vector<int> V_F1;
	vector<int> V_F2;
	vector<int> is_active;
	vector<int> AV_ID;
	vector<int>	boundary_vertexID;

	int		cellx_num;
	int		celly_num;
	double	x_min;
	double	x_max;
	double	y_min;
	double	y_max;
	double	lengthgrid_x;
	double	lengthgrid_y;
	vector<vector<int>>	cell_points;
};

