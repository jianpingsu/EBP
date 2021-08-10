#include "ShellData.h"
#include <iostream>

using namespace std;
using namespace Eigen;

void ShellData::update_shell()
{
	mv_num = m_V.rows();
	mf_num = m_T.rows();

	v_num = w_uv.rows();
	sf_num = s_T.rows();

	sv_num = v_num - mv_num;
	f_num = sf_num + mf_num;

	s_M = Eigen::VectorXd::Constant(sf_num, shell_factor);
}

void ShellData::mesh_improve()
{
	MatrixXd m_uv = w_uv.topRows(mv_num);
	MatrixXd V_bnd;
	V_bnd.resize(internal_bnd.size(), 2);
	for (int i = 0; i < internal_bnd.size(); i++)
	{
		V_bnd.row(i) = m_uv.row(internal_bnd(i));
	}

	if (frame_V.size() == 0)
	{
		VectorXd uv_max = m_uv.colwise().maxCoeff();
		VectorXd uv_min = m_uv.colwise().minCoeff();
		VectorXd uv_mid = (uv_max + uv_min) / 2.;
		double raduis1 = (uv_max(0) - uv_min(0)) / 2.;
		double raduis2 = 1.05 * raduis1;

		int frame_points = internal_bnd.size() / 10;
		if (frame_points > 200)
		{
			frame_points = 200;
		}
		frame_points = max(11, frame_points);

		double delta_angle = 2 * M_PI / frame_points;
		frame_V.resize(frame_points, 2);
		for (int i = 0; i < frame_points; ++i)
		{
			frame_V.row(i) << raduis2 * cos(i * delta_angle), raduis2 * sin(-i * delta_angle);
		}

		frame_ids = Eigen::VectorXi::LinSpaced(frame_V.rows(), mv_num, mv_num + frame_V.rows() - 1);
	}
	else
	{
		for (int i = 0; i < frame_V.rows(); ++i)
		{
			frame_V.row(i) = w_uv.row(frame_ids(i));
		}
	}

	MatrixXd V;
	MatrixXi E;

	V.resize(V_bnd.rows() + frame_V.rows(), V_bnd.cols());
	V << V_bnd, frame_V;

	E.resize(V.rows(), 2);
	for (int i = 0; i < E.rows(); i++)
		E.row(i) << i, i + 1;
	int acc_bs = 0;
	for (auto bs : bnd_sizes)
	{
		E(acc_bs + bs - 1, 1) = acc_bs;
		acc_bs += bs;
	}
	E(V.rows() - 1, 1) = acc_bs;
	assert(acc_bs == internal_bnd.size());

	//ycy H 存储的是每个chart的第一个三角形的重心
	MatrixXd H = MatrixXd::Zero(component_sizes.size(), 2);
	{
		int hole_f = 0;
		int hole_i = 0;
		for (auto cs : component_sizes) {
			for (int i = 0; i < 3; i++)
				H.row(hole_i) += m_uv.row(m_T(hole_f, i)); // redoing step 2
			hole_f += cs;
			hole_i++;
		}
	}
	H /= 3.;

	MatrixXd uv2;
	triangulate(V, E, H, uv2, s_T);
	auto bnd_n = internal_bnd.size();

	for (auto i = 0; i < s_T.rows(); i++)
	{
		for (auto j = 0; j < s_T.cols(); j++)
		{
			auto &x = s_T(i, j);
			if (x < bnd_n) x = internal_bnd(x);
			else x += m_uv.rows() - bnd_n;
		}
	}

	surface_F.resize(m_T.rows() + s_T.rows(), 3);
	surface_F << m_T, s_T;

	w_uv.conservativeResize(m_uv.rows() - bnd_n + uv2.rows(), 2);
	w_uv.bottomRows(uv2.rows() - bnd_n) = uv2.bottomRows(-bnd_n + uv2.rows());

	update_shell();
}

void ShellData::add_new_patch(const Eigen::MatrixXd &V_in, const Eigen::MatrixXi &F_ref, const Eigen::RowVectorXd &center)
{
	using namespace std;
	using namespace Eigen;

	VectorXd M;
	M.resize(F_ref.rows());
	{
		Eigen::Vector3d v0, v1, v2, e01, e02, normal_f;
		double area_f;
		for (size_t i = 0; i < F_ref.rows(); i++)
		{
			v0 = V_in.row(F_ref(i, 0));
			v1 = V_in.row(F_ref(i, 1));
			v2 = V_in.row(F_ref(i, 2));
			e01 = v1 - v0;
			e02 = v2 - v0;
			normal_f = e01.cross(e02);
			area_f = normal_f.norm() / 2.0;
			M(i) = area_f;
		}
	}

	m_M.conservativeResize(mf_num + M.size());
	m_M.bottomRows(M.size()) = M;
	mesh_measure += M.sum();

	const Eigen::MatrixXd& V_ref = V_in;
	Eigen::MatrixXd uv_init;
	Eigen::VectorXi bnd;
	Eigen::MatrixXd bnd_uv;
	std::vector<std::vector<int>> all_bnds;
	boundary_loop(F_ref, all_bnds);
	int num_holes = all_bnds.size() - 1;

	std::sort(all_bnds.begin(), all_bnds.end(), [](auto& a, auto&b)
	{
		return a.size() > b.size();
	});

	bnd = Map<Eigen::VectorXi>(all_bnds[0].data(), all_bnds[0].size());

	map_vertices_to_circle(V_ref, bnd, bnd_uv);
	bnd_uv *= sqrt(M.sum() / M_PI);
	bnd_uv.rowwise() += center;

	//std::cout << "Mesh Measure " << M.sum() << "; number holes " << num_holes << std::endl;

	if (num_holes == 0)
	{
		if (bnd.rows() == V_ref.rows())
		{
			std::cout << "All vert on boundary" << std::endl;
			uv_init.resize(V_ref.rows(), 2);
			for (int i = 0; i < bnd.rows(); i++)
			{
				uv_init.row(bnd(i)) = bnd_uv.row(i);
			}
		}
		else
		{
			Tutte(V_ref.rows(), F_ref, bnd, bnd_uv, uv_init);
		}
	}
	else
	{
		auto &F = F_ref;
		auto &V = V_in;
		auto &primary_bnd = bnd;
		// fill holes
		int n_filled_faces = 0;
		int real_F_num = F.rows();
		for (int i = 0; i < num_holes; i++)
		{
			n_filled_faces += all_bnds[i + 1].size();
		}
		MatrixXi F_filled(n_filled_faces + real_F_num, 3);
		F_filled.topRows(real_F_num) = F;

		int new_vert_id = V.rows();
		int new_face_id = real_F_num;

		for (int i = 0; i < num_holes; i++)
		{
			int cur_bnd_size = all_bnds[i + 1].size();
			auto it = all_bnds[i + 1].begin();
			auto back = all_bnds[i + 1].end() - 1;
			F_filled.row(new_face_id++) << *it, *back, new_vert_id;
			while (it != back)
			{
				F_filled.row(new_face_id++)
					<< *(it + 1), *(it), new_vert_id;
				it++;
			}
			new_vert_id++;
		}
		assert(new_face_id == F_filled.rows());
		assert(new_vert_id == V.rows() + num_holes);

		Tutte(V_ref.rows() + num_holes, F_filled, primary_bnd, bnd_uv, uv_init);
		uv_init.conservativeResize(V.rows(), 2);
	}

	//writeObj(uv_init, F_ref, "tutte.obj");

	component_sizes.push_back(F_ref.rows());

	if (mv_num == 0)
	{
		w_uv = uv_init;
	}
	else
	{
		MatrixXd m_uv = w_uv.topRows(mv_num);
		w_uv.resize(m_uv.rows() + uv_init.rows(), 2);
		w_uv << m_uv, uv_init;
	}

	for (auto cur_bnd : all_bnds)
	{
		internal_bnd.conservativeResize(internal_bnd.size() + cur_bnd.size());
		internal_bnd.bottomRows(cur_bnd.size()) = Map<ArrayXi>(cur_bnd.data(), cur_bnd.size()) + mv_num;
		bnd_sizes.push_back(cur_bnd.size());
	}

	m_T.conservativeResize(mf_num + F_ref.rows(), 3);
	m_T.bottomRows(F_ref.rows()) = F_ref.array() + mv_num;
	mf_num += F_ref.rows();

	m_V.conservativeResize(mv_num + V_ref.rows(), 3);
	m_V.bottomRows(V_ref.rows()) = V_ref;
	mv_num += V_ref.rows();

	frame_V.resize(0, 0);

	mesh_improve();
}

