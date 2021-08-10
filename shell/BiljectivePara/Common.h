#pragma once
#include"PardisoSolver.h"
#include<list>
#include <vector>
#include<set>
#include<fstream>
#include<math.h>
#include<Eigen/Core>

using namespace std;
void orderTheBoundary(vector<list<int>>& order_boundary, const vector<int>& nextlocation);
bool growFromP(int p, set<int>& isused, list<int>& order_boundary, const vector<int>& nextlocation);

void map_vertices_to_circle(const Eigen::MatrixXd& V, const Eigen::VectorXi& bnd, Eigen::MatrixXd& UV);

void Tutte(int V_N, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd, const Eigen::MatrixXd& bnd_uv, Eigen::MatrixXd & uv_init);
void preCalc_pardiso(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, PardisoSolver & pardiso);
void writeObj(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const string& outfile);
void boundary_loop(const Eigen::MatrixXi &F_ref, std::vector<std::vector<int>>& boundaryEdges);