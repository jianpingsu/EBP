#include "Parafun.h"

void Parafun::after_mesh_improve()
{
	total_num = d_.v_num;
	F_N = d_.f_num;
	V_N = d_.v_num;
	F0.resize(F_N);
	F1.resize(F_N);
	F2.resize(F_N);

	position_of_mesh.resize(2 * total_num);

	for (size_t i = 0; i < dim; i++)
	{
		position_of_mesh.block(i*total_num, 0, total_num, 1) = d_.w_uv.col(i);
	}

	init();
}
void Parafun::init()
{
	for (int i = 0; i < F_N; ++i)
	{
		F0[i] = d_.surface_F(i, 0);
		F1[i] = d_.surface_F(i, 1);
		F2[i] = d_.surface_F(i, 2);
	}

	handle_mintri();
	init_area();
	setvirtualtri();
	Pre_calculate();

	const double *pos = position_of_mesh.data();

	x_min = pos[0];
	x_max = pos[0];
	y_min = pos[V_N];
	y_max = pos[V_N];
	for (int i = 1; i < V_N; ++i)
	{
		if (pos[i] < x_min)
		{
			x_min = pos[i];
		}
		else if (pos[i] > x_max)
		{
			x_max = pos[i];
		}

		if (pos[i + V_N] < y_min)
		{
			y_min = pos[i + V_N];
		}
		else if (pos[i + V_N] > y_max)
		{
			y_max = pos[i + V_N];
		}
	}

	int frame_num = d_.frame_ids.size();
	average_length = 0;
	for (int j = 0; j < frame_num - 1; ++j)
	{
		average_length += (d_.frame_V.row(j) - d_.frame_V.row(j + 1)).norm();
	}
	average_length += (d_.frame_V.row(frame_num - 1) - d_.frame_V.row(0)).norm();
	average_length = average_length / frame_num;
	threhold = average_length / 8.0;
	threhold = 2 * threhold;	//dis * 2

	cellx_num = std::ceil((x_max - x_min) / average_length);
	celly_num = std::ceil((y_max - y_min) / average_length);
	cell_points.clear();
	cell_points.resize(cellx_num * celly_num);

	if (pardiso != NULL)
	{
		delete pardiso;
		pardiso = NULL;
	}
	pardiso = new PardisoSolver();
	pardiso->ia = pardiso_ia;
	pardiso->ja = pardiso_ja;
	pardiso->a.resize(pardiso_ja.size());
	pardiso->nnz = pardiso_ja.size();
	pardiso->num = 2 * V_N;
	pardiso->pardiso_init();

	fungrid(position_of_mesh);
	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double dis, E_b, energy2 = 0;
	for (int i = 0; i < AV_F_N; ++i)
	{
		f0 = V_F0[AV_ID[i]];
		f1 = V_F1[AV_ID[i]];
		f2 = V_F2[AV_ID[i]];

		x0 = pos[f0];	y0 = pos[f0 + V_N];

		x1 = pos[f1];	y1 = pos[f1 + V_N];

		x2 = pos[f2];	y2 = pos[f2 + V_N];

		dis = distance(x0, y0, x1, y1, x2, y2);
		if (dis < 0)
		{
			std::cout << "distance is zero" << std::endl;
		}
		E_b = (1 - threhold / dis) * (1 - threhold / dis);
		energy2 += E_b;
	}
	energy_barrier = energy2;

	if (is_first)
		is_first = false;
}
void Parafun::init_area()
{
	area.resize(F_N);
	int src_t_num = d_.m_T.rows();
	area_src.resize(src_t_num);
	for (int i = 0; i < src_t_num; ++i)
	{
		area_src[i] = d_.m_M(i);
		area[i] = d_.m_M(i);
	}
	for (int i = src_t_num; i < F_N; ++i)
	{
		area[i] = d_.s_M(i - src_t_num);
	}
}
void Parafun::setvirtualtri()

{
	VectorXi boundary_vertex = d_.frame_ids;
	BE_N = boundary_vertex.size();

	V_F_N = BE_N * (BE_N - 2);
	V_F0.resize(V_F_N);
	V_F1.resize(V_F_N);
	V_F2.resize(V_F_N);
	is_active.resize(V_F_N, -1);
	AV_ID.reserve(V_F_N);

	boundary_vertexID.clear();
	boundary_vertexID.resize(V_N, -1);

	int id = 0;
	for (int i = 0; i < BE_N; ++i)
	{
		boundary_vertexID[boundary_vertex(i)] = id;
		id++;
	}

	int id_start, id_end, id_mid;
	int k;
	for (int i = 0; i < BE_N; ++i)
	{
		id_start = boundary_vertex(i);
		id_end = boundary_vertex((i + 1) % BE_N);
		k = 0;

		for (int j = 0; j < BE_N; ++j)
		{
			if (id_start != boundary_vertex(j) && id_end != boundary_vertex(j))
			{
				V_F0[i * (BE_N - 2) + k] = id_start;
				V_F1[i * (BE_N - 2) + k] = id_end;
				V_F2[i * (BE_N - 2) + k] = boundary_vertex(j);
				++k;
			}
		}
	}
}
void Parafun::Pre_calculate()
{
	source_p00.resize(F_N);
	source_p01.resize(F_N);
	source_p10.resize(F_N);
	source_p11.resize(F_N);
	if (is_first)
	{
		for (int i = 0; i < d_.m_T.rows(); ++i)
		{
			double p00, p01, p10, p11;
			local_coordinate_inverse(i, p00, p01, p10, p11);
			source_p00[i] = p00;
			source_p01[i] = p01;
			source_p10[i] = p10;
			source_p11[i] = p11;
		}
	}

	for (int i = d_.m_T.rows(); i < F_N; ++i)
	{
		double p00, p01, p10, p11;
		local_coordinate_inverse_scaf(i, p00, p01, p10, p11);
		source_p00[i] = p00;
		source_p01[i] = p01;
		source_p10[i] = p10;
		source_p11[i] = p11;
	}
	update_p00 = source_p00;
	update_p01 = source_p01;
	update_p10 = source_p10;
	update_p11 = source_p11;

	pardiso_ia.clear(); pardiso_ia.reserve(2 * V_N + 1);
	pardiso_ja.clear(); pardiso_ja.reserve(8 * V_N + BE_N * BE_N * 2);
	typedef Triplet<int> T;
	std::vector<T> tripletlist;
	std::vector<std::set<int>> VV_tmp;
	VV_tmp.resize(V_N);
	for (size_t i = 0; i < d_.m_T.rows(); i++)
	{
		int vid[3];
		for (size_t j = 0; j < d_.m_T.cols(); j++)
		{
			vid[j] = d_.m_T(i, j);
		}
		VV_tmp[vid[0]].insert(vid[1]);
		VV_tmp[vid[0]].insert(vid[2]);
		VV_tmp[vid[1]].insert(vid[0]);
		VV_tmp[vid[1]].insert(vid[2]);
		VV_tmp[vid[2]].insert(vid[0]);
		VV_tmp[vid[2]].insert(vid[1]);
	}

	vector<int> s_vid;
	for (size_t i = 0; i < d_.s_T.rows(); i++)
	{
		s_vid.clear();
		for (size_t j = 0; j < d_.s_T.cols(); j++)
		{
			int s_id = d_.s_T(i, j);
			s_vid.push_back(s_id);
		}
		VV_tmp[s_vid[0]].insert(s_vid[1]);
		VV_tmp[s_vid[0]].insert(s_vid[2]);
		VV_tmp[s_vid[1]].insert(s_vid[0]);
		VV_tmp[s_vid[1]].insert(s_vid[2]);
		VV_tmp[s_vid[2]].insert(s_vid[0]);
		VV_tmp[s_vid[2]].insert(s_vid[1]);
	}

	vector<int> v_vid;
	for (int i = 0; i < d_.frame_ids.size(); ++i)
	{
		for (int j = 0; j < d_.frame_ids.size(); ++j)
		{
			if (j == i)
			{
				continue;
			}
			VV_tmp[d_.frame_ids(i)].insert(d_.frame_ids(j));
		}
	}

	for (int i = 0; i < V_N; i++)
	{
		pardiso_ia.push_back(pardiso_ja.size());
		VV_tmp[i].insert(i);
		vector<int> row_id;
		for (auto&var : VV_tmp[i])
		{
			row_id.push_back(var);
		}
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i);
		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			pardiso_ja.push_back(row_id[k]);
			tripletlist.push_back(T(i, row_id[k], dd));
			++dd;
		}
		for (int k = 0; k < row_id.size(); k++)
		{
			pardiso_ja.push_back(row_id[k] + V_N);
			tripletlist.push_back(T(i, row_id[k] + V_N, dd));
			++dd;
		}
	}
	for (int i = V_N; i < 2 * V_N; i++)
	{
		pardiso_ia.push_back(pardiso_ja.size());
		vector<int> row_id;
		for (auto&var : VV_tmp[i - V_N])
		{
			row_id.push_back(var);
		}
		vector<int>::iterator iter = std::find(row_id.begin(), row_id.end(), i - V_N);
		int dd = 0;
		for (int k = std::distance(row_id.begin(), iter); k < row_id.size(); k++)
		{
			pardiso_ja.push_back(row_id[k] + V_N);
			tripletlist.push_back(T(i, row_id[k] + V_N, dd));
			++dd;
		}
	}

	SparseMatrix<int> find_id_in_rows;
	find_id_in_rows.resize(2 * V_N, 2 * V_N);
	find_id_in_rows.setFromTriplets(tripletlist.begin(), tripletlist.end());
	pardiso_ia.push_back(pardiso_ja.size());
	id_h00.resize(F_N + V_F_N, -1); id_h01.resize(F_N + V_F_N, -1); id_h02.resize(F_N + V_F_N, -1); id_h03.resize(F_N + V_F_N, -1); id_h04.resize(F_N + V_F_N, -1); id_h05.resize(F_N + V_F_N, -1);
	id_h11.resize(F_N + V_F_N, -1); id_h12.resize(F_N + V_F_N, -1); id_h13.resize(F_N + V_F_N, -1); id_h14.resize(F_N + V_F_N, -1); id_h15.resize(F_N + V_F_N, -1);
	id_h22.resize(F_N + V_F_N, -1); id_h23.resize(F_N + V_F_N, -1); id_h24.resize(F_N + V_F_N, -1); id_h25.resize(F_N + V_F_N, -1);
	id_h33.resize(F_N + V_F_N, -1); id_h34.resize(F_N + V_F_N, -1); id_h35.resize(F_N + V_F_N, -1);
	id_h44.resize(F_N + V_F_N, -1); id_h45.resize(F_N + V_F_N, -1);
	id_h55.resize(F_N + V_F_N, -1);

	for (int i = 0; i < F_N; i++)
	{
		int f0 = F0[i]; int f1 = F1[i]; int f2 = F2[i]; int f3 = f0 + V_N; int f4 = f1 + V_N; int f5 = f2 + V_N;
		int min01 = min(f0, f1); int max01 = f0 + f1 - min01;
		int min02 = min(f0, f2); int max02 = f0 + f2 - min02;
		int min12 = min(f1, f2); int max12 = f1 + f2 - min12;
		id_h00[i] = pardiso_ia[f0]; id_h01[i] = pardiso_ia[min01] + find_id_in_rows.coeff(min01, max01); id_h02[i] = pardiso_ia[min02] + find_id_in_rows.coeff(min02, max02);
		id_h03[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f3); id_h04[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f4); id_h05[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f5);
		id_h11[i] = pardiso_ia[f1]; id_h12[i] = pardiso_ia[min12] + find_id_in_rows.coeff(min12, max12);
		id_h13[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f3); id_h14[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f4); id_h15[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f5);
		id_h22[i] = pardiso_ia[f2];
		id_h23[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f3); id_h24[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f4); id_h25[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f5);
		id_h33[i] = pardiso_ia[f3]; id_h34[i] = pardiso_ia[min01 + V_N] + find_id_in_rows.coeff(min01 + V_N, max01 + V_N); id_h35[i] = pardiso_ia[min02 + V_N] + find_id_in_rows.coeff(min02 + V_N, max02 + V_N);
		id_h44[i] = pardiso_ia[f4]; id_h45[i] = pardiso_ia[min12 + V_N] + find_id_in_rows.coeff(min12 + V_N, max12 + V_N);
		id_h55[i] = pardiso_ia[f5];

	}

	for (int i = F_N; i < F_N + V_F_N; i++)
	{
		int f0 = V_F0[i - F_N]; int f1 = V_F1[i - F_N]; int f2 = V_F2[i - F_N]; int f3 = V_F0[i - F_N] + V_N; int f4 = V_F1[i - F_N] + V_N; int f5 = V_F2[i - F_N] + V_N;
		int min01 = min(f0, f1); int max01 = f0 + f1 - min01;
		int min02 = min(f0, f2); int max02 = f0 + f2 - min02;
		int min12 = min(f1, f2); int max12 = f1 + f2 - min12;
		id_h00[i] = pardiso_ia[f0]; id_h01[i] = pardiso_ia[min01] + find_id_in_rows.coeff(min01, max01); id_h02[i] = pardiso_ia[min02] + find_id_in_rows.coeff(min02, max02);
		id_h03[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f3); id_h04[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f4); id_h05[i] = pardiso_ia[f0] + find_id_in_rows.coeff(f0, f5);
		id_h11[i] = pardiso_ia[f1]; id_h12[i] = pardiso_ia[min12] + find_id_in_rows.coeff(min12, max12);
		id_h13[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f3); id_h14[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f4); id_h15[i] = pardiso_ia[f1] + find_id_in_rows.coeff(f1, f5);
		id_h22[i] = pardiso_ia[f2];
		id_h23[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f3); id_h24[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f4); id_h25[i] = pardiso_ia[f2] + find_id_in_rows.coeff(f2, f5);
		id_h33[i] = pardiso_ia[f3]; id_h34[i] = pardiso_ia[min01 + V_N] + find_id_in_rows.coeff(min01 + V_N, max01 + V_N); id_h35[i] = pardiso_ia[min02 + V_N] + find_id_in_rows.coeff(min02 + V_N, max02 + V_N);
		id_h44[i] = pardiso_ia[f4]; id_h45[i] = pardiso_ia[min12 + V_N] + find_id_in_rows.coeff(min12 + V_N, max12 + V_N);
		id_h55[i] = pardiso_ia[f5];
	}
}
void Parafun::handle_mintri()
{
	double min_bnd_edge_len = numeric_limits<double>::infinity();
	int acc_bnd = 0;
	for (int i = 0; i < d_.bnd_sizes.size(); i++)
	{
		int current_size = d_.bnd_sizes[i];

		for (int e = acc_bnd; e < acc_bnd + current_size - 1; e++)
		{
			min_bnd_edge_len = (std::min)(min_bnd_edge_len, (d_.w_uv.row(d_.internal_bnd(e)) - d_.w_uv.row(d_.internal_bnd(e + 1))).squaredNorm());
		}
		min_bnd_edge_len = (std::min)(min_bnd_edge_len, (d_.w_uv.row(d_.internal_bnd(acc_bnd)) - d_.w_uv.row(d_.internal_bnd(acc_bnd + current_size - 1))).squaredNorm());
		acc_bnd += current_size;
	}

	area_threshold = min_bnd_edge_len / 4.0;
}

double Parafun::BPE(bool is_ip_convrate, bool is_slim_convrate)
{
	if (is_ip_convrate)
	{
		Update_source_same_t();
	}
	bool is_interp = is_ip_convrate && (Intp_T_Min < 0.999);
	bool is_slim = is_interp&& is_slim_convrate && (changetocm_flag < 0.99);

	if (is_slim)
	{
		SLIM(is_interp);
	}
	else
	{
		CM(is_interp);
	}
	d_.w_uv = Map<Matrix<double, -1, -1, Eigen::ColMajor>>(position_of_mesh.data(), total_num, dim);
	
	energy_all = energy_mesh + area.back() * energy_shell + barrer_coef * energy_barrier;
	return energy_all;
}
void Parafun::Update_source_same_t()
{
	double t_min = 1;
	int geqK = 0;
	int update_fn = d_.m_T.rows();
	vector<double> all_s0; all_s0.resize(update_fn);
	vector<double> all_s1; all_s1.resize(update_fn);
	vector<double> all_w00; all_w00.resize(update_fn);
	vector<double> all_w01; all_w01.resize(update_fn);
	vector<double> all_w10; all_w10.resize(update_fn);
	vector<double> all_w11; all_w11.resize(update_fn);

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det;
	double E_d;
	double tt;
	double new_sig0, new_sig1;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;
	double *position = position_of_mesh.data();

	for (int i = 0; i < update_fn; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];
		x0 = position[f0];
		y0 = position[f0 + total_num];
		x1 = position[f1];
		y1 = position[f1 + total_num];
		x2 = position[f2];
		y2 = position[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;
		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];
		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		det = j00 * j11 - j01 * j10;
		E_d = (1 + 1 / (det*det)) * (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);

		double alpha_0 = j00 + j11; double alpha_1 = j10 - j01;
		double beta_0 = j00 - j11; double beta_1 = j10 + j01;
		double alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1 * alpha_1);
		double beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1 * beta_1);

		double sig0 = alpha_norm + beta_norm;
		double sig1 = alpha_norm - beta_norm;
		all_s0[i] = sig0;
		all_s1[i] = sig1;

		if (beta_norm < 1e-15)
		{
			all_w00[i] = 0.0;
			all_w01[i] = 0.0;
			all_w10[i] = 0.0;
			all_w11[i] = 0.0;
		}
		else
		{
			double temp = 1 / (sig1*sig1 - sig0 * sig0);
			all_w00[i] = temp * (j00*j00 + j10 * j10 - 0.5*(sig0*sig0 + sig1 * sig1));
			all_w01[i] = temp * (j00*j01 + j10 * j11);
			all_w10[i] = temp * (j01*j00 + j11 * j10);
			all_w11[i] = temp * (j01*j01 + j11 * j11 - 0.5*(sig0*sig0 + sig1 * sig1));
		}

		if (E_d <= bound_distortion_K)
		{
			geqK++;
		}
		else
		{
			tt = newton_equation(sig0, sig1, bound_distortion_K);
			if (tt < t_min)
			{
				t_min = tt;
			}
		}
	}

	changetocm_flag = (double)geqK / update_fn;
	update_p00 = source_p00;
	update_p01 = source_p01;
	update_p10 = source_p10;
	update_p11 = source_p11;

	for (int i = 0; i < update_fn; ++i)
	{
		double sig0 = all_s0[i];
		double sig1 = all_s1[i];
		new_sig0 = pow(sig0, t_min - 1);
		new_sig1 = pow(sig1, t_min - 1);

		double delta_new = new_sig1 - new_sig0;
		double plus_new = 0.5*(new_sig1 + new_sig0);
		double w00 = delta_new * all_w00[i] + plus_new;
		double w01 = delta_new * all_w01[i];
		double w10 = delta_new * all_w10[i];
		double w11 = delta_new * all_w11[i] + plus_new;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];
		update_p00[i] = p00 * w00 + p01 * w10;
		update_p01[i] = p00 * w01 + p01 * w11;
		update_p10[i] = p10 * w00 + p11 * w10;
		update_p11[i] = p10 * w01 + p11 * w11;
	}

	Intp_T_Min = t_min;
}
void Parafun::SLIM(bool is_interp)
{
	double area_now;
	int f0, f1, f2;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	double x0, y0, x1, y1, x2, y2;

	double alpha_norm, beta_norm;

	double alpha_0, alpha_1, beta_0, beta_1;

	double sig0, sig1;

	double det, tr;
	double r0, r1, r2, r3;
	double d00, d01, d02,
		d10, d11, d12;

	double new_sig0, new_sig1;
	double temp;
	double w00, w01, w10, w11;
	double p1, p2, p3, w1, w2, w3;

	double h00, h01, h02, h03, h04, h05,
		h11, h12, h13, h14, h15,
		h22, h23, h24, h25,
		h33, h34, h35,
		h44, h45,
		h55;
	double *position = position_of_mesh.data();

	int nnz = pardiso_ja.size();
	pardiso_a.clear(); pardiso_b.clear();
	pardiso_a.resize(nnz, 0.0);
	pardiso_b.resize(2 * V_N, 0.0);

	double* tmp_p00;
	double* tmp_p01;
	double* tmp_p10;
	double* tmp_p11;

	if (is_interp)
	{
		tmp_p00 = update_p00.data();
		tmp_p01 = update_p01.data();
		tmp_p10 = update_p10.data();
		tmp_p11 = update_p11.data();
	}
	else
	{
		tmp_p00 = source_p00.data();
		tmp_p01 = source_p01.data();
		tmp_p10 = source_p10.data();
		tmp_p11 = source_p11.data();
	}

	int src_t_num = d_.m_T.rows();
	for (int i = 0; i < src_t_num; ++i)
	{
		area_now = area[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];
		x0 = position[f0];
		y0 = position[f0 + total_num];
		x1 = position[f1];
		y1 = position[f1 + total_num];
		x2 = position[f2];
		y2 = position[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;
		p00 = tmp_p00[i]; p01 = tmp_p01[i]; p10 = tmp_p10[i]; p11 = tmp_p11[i];
		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11; beta_1 = j10 + j01;
		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1*alpha_1);
		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1*beta_1);

		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;
		new_sig0 = sqrt(1 + 1 / sig0 + 1 / (sig0*sig0) + 1 / (sig0*sig0*sig0)); new_sig1 = sqrt(1 + 1 / sig1 + 1 / (sig1*sig1) + 1 / (sig1*sig1*sig1));

		if (abs(sig1 - sig0) < 1e-10)
		{
			temp = 0;
		}
		else
		{
			temp = (new_sig1 - new_sig0) / (sig1*sig1 - sig0 * sig0);
		}

		w00 = temp*(j00*j00 + j01*j01 - 0.5*(sig0*sig0 + sig1*sig1)) + 0.5*(new_sig0 + new_sig1);
		w01 = temp*(j00*j10 + j01*j11);
		w10 = temp*(j10*j00 + j11*j01);
		w11 = temp*(j10*j10 + j11*j11 - 0.5*(sig0*sig0 + sig1*sig1)) + 0.5*(new_sig0 + new_sig1);
		p1 = p00*p00 + p01*p01; p2 = p00*p10 + p01*p11; p3 = p10*p10 + p11*p11;
		w1 = w00*w00 + w10*w10; w2 = w00*w01 + w10*w11; w3 = w01*w01 + w11*w11;

		h00 = area_now *(p1 + p2 + p2 + p3)*w1; h01 = -area_now *(p1 + p2)*w1; h02 = -area_now *(p2 + p3)*w1; h03 = area_now *(p1 + p2 + p2 + p3)*w2; h04 = -area_now *(p1 + p2)*w2; h05 = -area_now *(p2 + p3)*w2;
		h11 = area_now *p1*w1;                  h12 = area_now *p2*w1;    	 h13 = -area_now *(p1 + p2)*w2; h14 = area_now *p1*w2;                  h15 = area_now *p2*w2;
		h22 = area_now *p3*w1;                  h23 = -area_now *(p2 + p3)*w2; h24 = area_now *p2*w2;         h25 = area_now *p3*w2;
		h33 = area_now *(p1 + p2 + p2 + p3)*w3; h34 = -area_now *(p1 + p2)*w3; h35 = -area_now *(p2 + p3)*w3;
		h44 = area_now *p1*w3;                  h45 = area_now *p2*w3;
		h55 = area_now *p3*w3;

		det = j00*j11 - j01*j10;
		tr = (j00*j00 + j01*j01 + j10*j10 + j11*j11);
		d00 = -p00 - p10; d01 = p00; d02 = p10;
		d10 = -p01 - p11; d11 = p01; d12 = p11;
		r0 = area_now * ((1 + 1 / (det*det))*j00 - tr*j11 / (det*det*det));
		r1 = area_now * ((1 + 1 / (det*det))*j01 + tr*j10 / (det*det*det));
		r2 = area_now * ((1 + 1 / (det*det))*j10 + tr*j01 / (det*det*det));
		r3 = area_now * ((1 + 1 / (det*det))*j11 - tr*j00 / (det*det*det));

		pardiso_b[f0] -= r0*d00 + r1*d10;
		pardiso_b[f1] -= r0*d01 + r1*d11;
		pardiso_b[f2] -= r0*d02 + r1*d12;
		pardiso_b[f0 + V_N] -= r2*d00 + r3*d10;
		pardiso_b[f1 + V_N] -= r2*d01 + r3*d11;
		pardiso_b[f2 + V_N] -= r2*d02 + r3*d12;

		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;
	}

	for (int i = src_t_num; i < F_N; ++i)
	{
		area_now = area[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];
		x0 = position[f0];
		y0 = position[f0 + total_num];
		x1 = position[f1];
		y1 = position[f1 + total_num];
		x2 = position[f2];
		y2 = position[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;
		p00 = tmp_p00[i]; p01 = tmp_p01[i]; p10 = tmp_p10[i]; p11 = tmp_p11[i];
		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;
		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11; beta_1 = j10 + j01;

		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1*alpha_1);
		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1*beta_1);
		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;
		new_sig0 = sqrt(1 + 1 / sig0 + 1 / (sig0*sig0) + 1 / (sig0*sig0*sig0)); new_sig1 = sqrt(1 + 1 / sig1 + 1 / (sig1*sig1) + 1 / (sig1*sig1*sig1));

		if (abs(sig1 - sig0) < 1e-10)
		{
			temp = 0;
		}
		else
		{
			temp = (new_sig1 - new_sig0) / (sig1*sig1 - sig0 * sig0);
		}

		w00 = temp*(j00*j00 + j01*j01 - 0.5*(sig0*sig0 + sig1*sig1)) + 0.5*(new_sig0 + new_sig1);
		w01 = temp*(j00*j10 + j01*j11);
		w10 = temp*(j10*j00 + j11*j01);
		w11 = temp*(j10*j10 + j11*j11 - 0.5*(sig0*sig0 + sig1*sig1)) + 0.5*(new_sig0 + new_sig1);
		p1 = p00*p00 + p01*p01; p2 = p00*p10 + p01*p11; p3 = p10*p10 + p11*p11;
		w1 = w00*w00 + w10*w10; w2 = w00*w01 + w10*w11; w3 = w01*w01 + w11*w11;

		h00 = area_now *(p1 + p2 + p2 + p3)*w1; h01 = -area_now *(p1 + p2)*w1; h02 = -area_now *(p2 + p3)*w1; h03 = area_now *(p1 + p2 + p2 + p3)*w2; h04 = -area_now *(p1 + p2)*w2; h05 = -area_now *(p2 + p3)*w2;
		h11 = area_now *p1*w1;                  h12 = area_now *p2*w1;    	 h13 = -area_now *(p1 + p2)*w2; h14 = area_now *p1*w2;                  h15 = area_now *p2*w2;
		h22 = area_now *p3*w1;                  h23 = -area_now *(p2 + p3)*w2; h24 = area_now *p2*w2;         h25 = area_now *p3*w2;
		h33 = area_now *(p1 + p2 + p2 + p3)*w3; h34 = -area_now *(p1 + p2)*w3; h35 = -area_now *(p2 + p3)*w3;
		h44 = area_now *p1*w3;                  h45 = area_now *p2*w3;
		h55 = area_now *p3*w3;

		det = j00*j11 - j01*j10;
		tr = (j00*j00 + j01*j01 + j10*j10 + j11*j11);
		d00 = -p00 - p10; d01 = p00; d02 = p10;
		d10 = -p01 - p11; d11 = p01; d12 = p11;
		r0 = area_now * ((1 + 1 / (det*det))*j00 - tr*j11 / (det*det*det));
		r1 = area_now * ((1 + 1 / (det*det))*j01 + tr*j10 / (det*det*det));
		r2 = area_now * ((1 + 1 / (det*det))*j10 + tr*j01 / (det*det*det));
		r3 = area_now * ((1 + 1 / (det*det))*j11 - tr*j00 / (det*det*det));

		pardiso_b[f0] -= r0*d00 + r1*d10;
		pardiso_b[f1] -= r0*d01 + r1*d11;
		pardiso_b[f2] -= r0*d02 + r1*d12;
		pardiso_b[f0 + V_N] -= r2*d00 + r3*d10;
		pardiso_b[f1 + V_N] -= r2*d01 + r3*d11;
		pardiso_b[f2 + V_N] -= r2*d02 + r3*d12;

		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;
	}

	double len01, len0i, len1i, coef1, coef2, dis;
	int i;
	double s0, s1, e0, e1, p0;
	double h_u, uu, aa, bb, cc;
	double a1x0, a1x1, a1x2, a1x3, a1x4, a1x5,
		a2x0, a2x1, a2x2, a2x3, a2x4, a2x5;
	for (int ii = 0; ii < AV_F_N; ++ii)
	{
		f0 = V_F0[AV_ID[ii]];
		f1 = V_F1[AV_ID[ii]];
		f2 = V_F2[AV_ID[ii]];
		x0 = position[f0];
		y0 = position[f0 + V_N];
		x1 = position[f1];
		y1 = position[f1 + V_N];
		x2 = position[f2];
		y2 = position[f2 + V_N];

		s0 = x0 + x1;	s1 = y0 + y1;
		e0 = x1 + x0;		e1 = y1 + y0;
		p0 = x2;	p1 = y2;
		len01 = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
		len0i = sqrt((p0 - s0)*(p0 - s0) + (p1 - s1)*(p1 - s1));
		len1i = sqrt((p0 - e0)*(p0 - e0) + (p1 - e1)*(p1 - e1));
		dis = len0i + len1i - len01;

		h_u = 2 * barrer_coef * threhold * (dis - threhold) / dis / dis / dis;
		uu = 2 * barrer_coef * threhold * (3 * threhold - 2 * dis) / dis / dis / dis / dis;

		a1x0 = (x0 - x2) / len0i - (x0 - x1) / len01;
		a1x1 = (x1 - x2) / len1i - (x1 - x0) / len01;
		a1x2 = (x2 - x0) / len0i + (x2 - x1) / len1i;
		a1x3 = (y0 - y2) / len0i - (y0 - y1) / len01;
		a1x4 = (y1 - y2) / len1i - (y1 - y0) / len01;
		a1x5 = (y2 - y0) / len0i + (y2 - y1) / len1i;

		coef1 = len01 * len01;
		coef2 = h_u / coef1 / len01;

		aa = (x0 - x1)*(x0 - x1);
		bb = (x1 - x0)*(y1 - y0);
		cc = (y0 - y1)*(y0 - y1);
		h00 = coef2 * (aa - coef1); h01 = coef2 * (coef1 - aa); h02 = 0.0; h03 = coef2 * bb; h04 = -coef2 * bb; h05 = 0.0;
		h11 = coef2 * (aa - coef1); h12 = 0.0; h13 = -coef2 * bb; h14 = coef2 * bb; h15 = 0.0;
		h22 = 0.0; h23 = 0.0; h24 = 0.0; h25 = 0.0;
		h33 = coef2 * (cc - coef1); h34 = coef2 * (coef1 - cc); h35 = 0.0;
		h44 = coef2 * (cc - coef1); h45 = 0.0;
		h55 = 0.0;

		h00 += uu*a1x0*a1x0; h01 += uu*a1x0*a1x1; h02 += uu*a1x0*a1x2; h03 += uu*a1x0*a1x3; h04 += uu*a1x0*a1x4; h05 += uu*a1x0*a1x5;
		h11 += uu*a1x1*a1x1; h12 += uu*a1x1*a1x2; h13 += uu*a1x1*a1x3; h14 += uu*a1x1*a1x4; h15 += uu*a1x1*a1x5;
		h22 += uu*a1x2*a1x2; h23 += uu*a1x2*a1x3; h24 += uu*a1x2*a1x4; h25 += uu*a1x2*a1x5;
		h33 += uu*a1x3*a1x3; h34 += uu*a1x3*a1x4; h35 += uu*a1x3*a1x5;
		h44 += uu*a1x4*a1x4; h45 += uu*a1x4*a1x5;
		h55 += uu*a1x5*a1x5;

		pardiso_b[f0] -= h_u * (a1x0);
		pardiso_b[f1] -= h_u * (a1x1);
		pardiso_b[f2] -= h_u * (a1x2);
		pardiso_b[f0 + V_N] -= h_u * (a1x3);
		pardiso_b[f1 + V_N] -= h_u * (a1x4);
		pardiso_b[f2 + V_N] -= h_u * (a1x5);

		i = F_N + AV_ID[ii];
		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;

	}

	pardiso_a[0] += 1.0;
	pardiso_a[pardiso_ia[V_N]] += 1.0;
	pardiso->a = pardiso_a;
	pardiso->rhs = pardiso_b;
	long time_beg, time_end;
	time_beg = clock();
	pardiso->factorize();
	time_end = clock();
	double time_consumption = (time_end - time_beg) / 1000.0;
	//std::cout << "factorize" << time_consumption << std::endl;
	pardiso->pardiso_solver();

	vector<double> result_d = pardiso->result;
	VectorXd negative_grad(2 * total_num), d(2 * total_num);
	for (int i = 0; i < V_N; i++)
	{
		negative_grad(i) = pardiso_b[i];
		negative_grad(i + total_num) = pardiso_b[i + V_N];
		d(i) = result_d[i];
		d(i + total_num) = result_d[i + V_N];
	}
	pardiso->free_numerical_factorization_memory();

	double temp_t;
	max_step(position_of_mesh, d, temp_t);
	double alpha = min(1.0, 0.8 * temp_t);
	backtracking_line_search(position_of_mesh, d, negative_grad, alpha, is_interp);
	position_of_mesh += alpha * d;
	Energysource();
}
void Parafun::CM(bool is_interp)
{
	double area_now;
	int f0, f1, f2;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	double x0, y0, x1, y1, x2, y2;
	double hi_0, hi_1;
	double alpha_0, alpha_1, beta_0, beta_1;
	double ss1, ss2, sig0, sig1;
	double alpha_norm, beta_norm;
	double h_u, h_v, walpha, wbeta;
	double a1x0, a1x1, a1x2, a1x3, a1x4, a1x5,
		a2x0, a2x1, a2x2, a2x3, a2x4, a2x5;

	double aa, bb, cc;
	double uu, vv, uv;
	double u, v;
	double h00, h01, h02, h03, h04, h05,
		h11, h12, h13, h14, h15,
		h22, h23, h24, h25,
		h33, h34, h35,
		h44, h45,
		h55;
	double *position = position_of_mesh.data();
	int nnz = pardiso_ja.size();
	pardiso_a.clear(); pardiso_b.clear();
	pardiso_a.resize(nnz, 0.0);
	pardiso_b.resize(2 * V_N, 0.0);

	double* tmp_p00;
	double* tmp_p01;
	double* tmp_p10;
	double* tmp_p11;
	if (is_interp)
	{
		tmp_p00 = update_p00.data();
		tmp_p01 = update_p01.data();
		tmp_p10 = update_p10.data();
		tmp_p11 = update_p11.data();
	}
	else
	{
		tmp_p00 = source_p00.data();
		tmp_p01 = source_p01.data();
		tmp_p10 = source_p10.data();
		tmp_p11 = source_p11.data();
	}

	for (int i = 0; i < d_.m_T.rows(); ++i)
	{
		area_now = area[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];
		x0 = position[f0];
		y0 = position[f0 + V_N];
		x1 = position[f1];
		y1 = position[f1 + V_N];
		x2 = position[f2];
		y2 = position[f2 + V_N];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;
		p00 = tmp_p00[i]; p01 = tmp_p01[i]; p10 = tmp_p10[i]; p11 = tmp_p11[i];
		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;

		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11;  beta_1 = j10 + j01;
		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1*alpha_1);
		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1*beta_1);
		ss1 = (p00)*(p00 + p10) + (p01)*(p01 + p11);
		ss2 = (p10)*(p00 + p10) + (p11)*(p01 + p11);

		double h1 = p00*p00 + p01*p01;
		double h2 = p00*p10 + p01*p11;
		double h3 = p10*p10 + p11*p11;
		double h4 = p00*p11 - p01*p10;

		a1x0 = alpha_0*(-p00 - p10) + alpha_1*(p01 + p11);  a1x1 = alpha_0*p00 - alpha_1*p01; a1x2 = alpha_0*p10 - alpha_1*p11;
		a1x3 = alpha_0*(-p01 - p11) + alpha_1*(-p00 - p10); a1x4 = alpha_0*p01 + alpha_1*p00; a1x5 = alpha_0*p11 + alpha_1*p10;
		a2x0 = beta_0*(-p00 - p10) + beta_1*(-p01 - p11);   a2x1 = beta_0*p00 + beta_1*p01;   a2x2 = beta_0*p10 + beta_1*p11;
		a2x3 = beta_0*(p01 + p11) + beta_1*(-p00 - p10);    a2x4 = -beta_0*p01 + beta_1*p00;  a2x5 = -beta_0*p11 + beta_1*p10;
		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		hi_0 = 2 + 6 * 1 / (sig0*sig0*sig0*sig0); hi_1 = 2 + 6 * 1 / (sig1*sig1*sig1*sig1);
		aa = 0.25 / alpha_norm; bb = 0.25 / beta_norm;
		uu = aa*aa*(area_now*hi_0 + area_now*hi_1);
		vv = bb*bb*(area_now*hi_0 + area_now*hi_1);
		uv = aa*bb*(area_now*hi_0 - area_now*hi_1);
		h_u = area_now * (2 * sig0 - 2 * 1 / (sig0*sig0*sig0));
		h_v = area_now * (2 * sig1 - 2 * 1 / (sig1*sig1*sig1));

		walpha = h_u + h_v;
		wbeta = h_u - h_v;
		double hwa1 = (walpha * 0.25 / alpha_norm); double hwa2 = -(walpha * 0.25*0.25 / (alpha_norm*alpha_norm*alpha_norm));
		double hwb1 = (wbeta * 0.25 / beta_norm); double hwb2 = -(wbeta *0.25*0.25 / (beta_norm*beta_norm*beta_norm));

		h00 = uu*a1x0*a1x0 + vv*a2x0*a2x0 + uv*a1x0*a2x0 + uv*a2x0*a1x0; h01 = uu*a1x0*a1x1 + vv*a2x0*a2x1 + uv*a1x0*a2x1 + uv*a2x0*a1x1; h02 = uu*a1x0*a1x2 + vv*a2x0*a2x2 + uv*a1x0*a2x2 + uv*a2x0*a1x2; h03 = uu*a1x0*a1x3 + vv*a2x0*a2x3 + uv*a1x0*a2x3 + uv*a2x0*a1x3; h04 = uu*a1x0*a1x4 + vv*a2x0*a2x4 + uv*a1x0*a2x4 + uv*a2x0*a1x4; h05 = uu*a1x0*a1x5 + vv*a2x0*a2x5 + uv*a1x0*a2x5 + uv*a2x0*a1x5;
		h11 = uu*a1x1*a1x1 + vv*a2x1*a2x1 + uv*a1x1*a2x1 + uv*a2x1*a1x1; h12 = uu*a1x1*a1x2 + vv*a2x1*a2x2 + uv*a1x1*a2x2 + uv*a2x1*a1x2; h13 = uu*a1x1*a1x3 + vv*a2x1*a2x3 + uv*a1x1*a2x3 + uv*a2x1*a1x3; h14 = uu*a1x1*a1x4 + vv*a2x1*a2x4 + uv*a1x1*a2x4 + uv*a2x1*a1x4; h15 = uu*a1x1*a1x5 + vv*a2x1*a2x5 + uv*a1x1*a2x5 + uv*a2x1*a1x5;
		h22 = uu*a1x2*a1x2 + vv*a2x2*a2x2 + uv*a1x2*a2x2 + uv*a2x2*a1x2; h23 = uu*a1x2*a1x3 + vv*a2x2*a2x3 + uv*a1x2*a2x3 + uv*a2x2*a1x3; h24 = uu*a1x2*a1x4 + vv*a2x2*a2x4 + uv*a1x2*a2x4 + uv*a2x2*a1x4; h25 = uu*a1x2*a1x5 + vv*a2x2*a2x5 + uv*a1x2*a2x5 + uv*a2x2*a1x5;
		h33 = uu*a1x3*a1x3 + vv*a2x3*a2x3 + uv*a1x3*a2x3 + uv*a2x3*a1x3; h34 = uu*a1x3*a1x4 + vv*a2x3*a2x4 + uv*a1x3*a2x4 + uv*a2x3*a1x4; h35 = uu*a1x3*a1x5 + vv*a2x3*a2x5 + uv*a1x3*a2x5 + uv*a2x3*a1x5;
		h44 = uu*a1x4*a1x4 + vv*a2x4*a2x4 + uv*a1x4*a2x4 + uv*a2x4*a1x4; h45 = uu*a1x4*a1x5 + vv*a2x4*a2x5 + uv*a1x4*a2x5 + uv*a2x4*a1x5;
		h55 = uu*a1x5*a1x5 + vv*a2x5*a2x5 + uv*a1x5*a2x5 + uv*a2x5*a1x5;

		if (walpha >= 0)
		{
			h00 += hwa1*(ss1 + ss2) + hwa2*a1x0*a1x0; h01 += hwa1*(-ss1) + hwa2*a1x0*a1x1; h02 += hwa1*(-ss2) + hwa2*a1x0*a1x2; h03 += hwa2*a1x0*a1x3; h04 += hwa1*(h4)+hwa2*a1x0*a1x4; h05 += hwa1*(-h4) + hwa2*a1x0*a1x5;
			h11 += hwa1*(h1)+hwa2*a1x1*a1x1;        h12 += hwa1*(h2)+hwa2*a1x1*a1x2;    h13 += hwa1*(-h4) + hwa2*a1x1*a1x3; h14 += hwa2*a1x1*a1x4; h15 += hwa1*(h4)+hwa2*a1x1*a1x5;
			h22 += hwa1*(h3)+hwa2*a1x2*a1x2;        h23 += hwa1*(h4)+hwa2*a1x2*a1x3;    h24 += hwa1*(-h4) + hwa2*a1x2*a1x4; h25 += hwa2*a1x2*a1x5;
			h33 += hwa1*(ss1 + ss2) + hwa2*a1x3*a1x3; h34 += hwa1*(-ss1) + hwa2*a1x3*a1x4; h35 += hwa1*(-ss2) + hwa2*a1x3*a1x5;
			h44 += hwa1*(h1)+hwa2*a1x4*a1x4;        h45 += hwa1*(h2)+hwa2*a1x4*a1x5;
			h55 += hwa1*(h3)+hwa2*a1x5*a1x5;

		}
		h00 += hwb1*(ss1 + ss2) + hwb2*a2x0*a2x0; h01 += hwb1*(-ss1) + hwb2*a2x0*a2x1; h02 += hwb1*(-ss2) + hwb2*a2x0*a2x2; h03 += hwb2*a2x0*a2x3; h04 += hwb1*(-h4) + hwb2*a2x0*a2x4; h05 += hwb1*(h4)+hwb2*a2x0*a2x5;
		h11 += hwb1*(h1)+hwb2*a2x1*a2x1;        h12 += hwb1*(h2)+hwb2*a2x1*a2x2;    h13 += hwb1*(h4)+hwb2*a2x1*a2x3;    h14 += hwb2*a2x1*a2x4; h15 += hwb1*(-h4) + hwb2*a2x1*a2x5;
		h22 += hwb1*(h3)+hwb2*a2x2*a2x2;        h23 += hwb1*(-h4) + hwb2*a2x2*a2x3; h24 += hwb1*(h4)+hwb2*a2x2*a2x4;    h25 += hwb2*a2x2*a2x5;
		h33 += hwb1*(ss1 + ss2) + hwb2*a2x3*a2x3; h34 += hwb1*(-ss1) + hwb2*a2x3*a2x4; h35 += hwb1*(-ss2) + hwb2*a2x3*a2x5;
		h44 += hwb1*(h1)+hwb2*a2x4*a2x4;        h45 += hwb1*(h2)+hwb2*a2x4*a2x5;
		h55 += hwb1*(h3)+hwb2*a2x5*a2x5;

		u = aa*walpha; v = bb*wbeta;
		pardiso_b[f0] -= (u*a1x0 + v*a2x0);
		pardiso_b[f1] -= (u*a1x1 + v*a2x1);
		pardiso_b[f2] -= (u*a1x2 + v*a2x2);
		pardiso_b[f0 + V_N] -= (u*a1x3 + v*a2x3);
		pardiso_b[f1 + V_N] -= (u*a1x4 + v*a2x4);
		pardiso_b[f2 + V_N] -= (u*a1x5 + v*a2x5);

		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;
	}

	for (int i = d_.m_T.rows(); i < F_N; i++)
	{
		area_now = area[i];
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];
		x0 = position[f0];
		y0 = position[f0 + V_N];
		x1 = position[f1];
		y1 = position[f1 + V_N];
		x2 = position[f2];
		y2 = position[f2 + V_N];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;
		p00 = tmp_p00[i]; p01 = tmp_p01[i]; p10 = tmp_p10[i]; p11 = tmp_p11[i];
		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;
		alpha_0 = j00 + j11; alpha_1 = j10 - j01;
		beta_0 = j00 - j11;  beta_1 = j10 + j01;

		alpha_norm = 0.5*sqrt(alpha_0*alpha_0 + alpha_1*alpha_1);
		beta_norm = 0.5*sqrt(beta_0*beta_0 + beta_1*beta_1);
		if (beta_norm < 1e-30)
		{
			beta_norm = 1e-10;
		}

		ss1 = (p00)*(p00 + p10) + (p01)*(p01 + p11);
		ss2 = (p10)*(p00 + p10) + (p11)*(p01 + p11);
		double h1 = p00*p00 + p01*p01;
		double h2 = p00*p10 + p01*p11;
		double h3 = p10*p10 + p11*p11;
		double h4 = p00*p11 - p01*p10;

		a1x0 = alpha_0*(-p00 - p10) + alpha_1*(p01 + p11);  a1x1 = alpha_0*p00 - alpha_1*p01; a1x2 = alpha_0*p10 - alpha_1*p11;
		a1x3 = alpha_0*(-p01 - p11) + alpha_1*(-p00 - p10); a1x4 = alpha_0*p01 + alpha_1*p00; a1x5 = alpha_0*p11 + alpha_1*p10;
		a2x0 = beta_0*(-p00 - p10) + beta_1*(-p01 - p11);   a2x1 = beta_0*p00 + beta_1*p01;   a2x2 = beta_0*p10 + beta_1*p11;
		a2x3 = beta_0*(p01 + p11) + beta_1*(-p00 - p10);    a2x4 = -beta_0*p01 + beta_1*p00;  a2x5 = -beta_0*p11 + beta_1*p10;
		sig0 = alpha_norm + beta_norm;
		sig1 = alpha_norm - beta_norm;

		hi_0 = 2 + 6 * 1 / (sig0*sig0*sig0*sig0); hi_1 = 2 + 6 * 1 / (sig1*sig1*sig1*sig1);
		aa = 0.25 / alpha_norm; bb = 0.25 / beta_norm;
		uu = aa*aa*(area_now*hi_0 + area_now*hi_1);
		vv = bb*bb*(area_now*hi_0 + area_now*hi_1);
		uv = aa*bb*(area_now*hi_0 - area_now*hi_1);
		h_u = area_now * (2 * sig0 - 2 * 1 / (sig0*sig0*sig0));
		h_v = area_now * (2 * sig1 - 2 * 1 / (sig1*sig1*sig1));
		walpha = h_u + h_v;
		wbeta = h_u - h_v;

		double hwa1 = (walpha * 0.25 / alpha_norm); double hwa2 = -(walpha * 0.25*0.25 / (alpha_norm*alpha_norm*alpha_norm));
		double hwb1 = (wbeta * 0.25 / beta_norm); double hwb2 = -(wbeta *0.25*0.25 / (beta_norm*beta_norm*beta_norm));
		h00 = uu*a1x0*a1x0 + vv*a2x0*a2x0 + uv*a1x0*a2x0 + uv*a2x0*a1x0; h01 = uu*a1x0*a1x1 + vv*a2x0*a2x1 + uv*a1x0*a2x1 + uv*a2x0*a1x1; h02 = uu*a1x0*a1x2 + vv*a2x0*a2x2 + uv*a1x0*a2x2 + uv*a2x0*a1x2; h03 = uu*a1x0*a1x3 + vv*a2x0*a2x3 + uv*a1x0*a2x3 + uv*a2x0*a1x3; h04 = uu*a1x0*a1x4 + vv*a2x0*a2x4 + uv*a1x0*a2x4 + uv*a2x0*a1x4; h05 = uu*a1x0*a1x5 + vv*a2x0*a2x5 + uv*a1x0*a2x5 + uv*a2x0*a1x5;
		h11 = uu*a1x1*a1x1 + vv*a2x1*a2x1 + uv*a1x1*a2x1 + uv*a2x1*a1x1; h12 = uu*a1x1*a1x2 + vv*a2x1*a2x2 + uv*a1x1*a2x2 + uv*a2x1*a1x2; h13 = uu*a1x1*a1x3 + vv*a2x1*a2x3 + uv*a1x1*a2x3 + uv*a2x1*a1x3; h14 = uu*a1x1*a1x4 + vv*a2x1*a2x4 + uv*a1x1*a2x4 + uv*a2x1*a1x4; h15 = uu*a1x1*a1x5 + vv*a2x1*a2x5 + uv*a1x1*a2x5 + uv*a2x1*a1x5;
		h22 = uu*a1x2*a1x2 + vv*a2x2*a2x2 + uv*a1x2*a2x2 + uv*a2x2*a1x2; h23 = uu*a1x2*a1x3 + vv*a2x2*a2x3 + uv*a1x2*a2x3 + uv*a2x2*a1x3; h24 = uu*a1x2*a1x4 + vv*a2x2*a2x4 + uv*a1x2*a2x4 + uv*a2x2*a1x4; h25 = uu*a1x2*a1x5 + vv*a2x2*a2x5 + uv*a1x2*a2x5 + uv*a2x2*a1x5;
		h33 = uu*a1x3*a1x3 + vv*a2x3*a2x3 + uv*a1x3*a2x3 + uv*a2x3*a1x3; h34 = uu*a1x3*a1x4 + vv*a2x3*a2x4 + uv*a1x3*a2x4 + uv*a2x3*a1x4; h35 = uu*a1x3*a1x5 + vv*a2x3*a2x5 + uv*a1x3*a2x5 + uv*a2x3*a1x5;
		h44 = uu*a1x4*a1x4 + vv*a2x4*a2x4 + uv*a1x4*a2x4 + uv*a2x4*a1x4; h45 = uu*a1x4*a1x5 + vv*a2x4*a2x5 + uv*a1x4*a2x5 + uv*a2x4*a1x5;
		h55 = uu*a1x5*a1x5 + vv*a2x5*a2x5 + uv*a1x5*a2x5 + uv*a2x5*a1x5;

		if (walpha >= 0)
		{
			h00 += hwa1*(ss1 + ss2) + hwa2*a1x0*a1x0; h01 += hwa1*(-ss1) + hwa2*a1x0*a1x1; h02 += hwa1*(-ss2) + hwa2*a1x0*a1x2; h03 += hwa2*a1x0*a1x3; h04 += hwa1*(h4)+hwa2*a1x0*a1x4; h05 += hwa1*(-h4) + hwa2*a1x0*a1x5;
			h11 += hwa1*(h1)+hwa2*a1x1*a1x1;        h12 += hwa1*(h2)+hwa2*a1x1*a1x2;    h13 += hwa1*(-h4) + hwa2*a1x1*a1x3; h14 += hwa2*a1x1*a1x4; h15 += hwa1*(h4)+hwa2*a1x1*a1x5;
			h22 += hwa1*(h3)+hwa2*a1x2*a1x2;        h23 += hwa1*(h4)+hwa2*a1x2*a1x3;    h24 += hwa1*(-h4) + hwa2*a1x2*a1x4; h25 += hwa2*a1x2*a1x5;
			h33 += hwa1*(ss1 + ss2) + hwa2*a1x3*a1x3; h34 += hwa1*(-ss1) + hwa2*a1x3*a1x4; h35 += hwa1*(-ss2) + hwa2*a1x3*a1x5;
			h44 += hwa1*(h1)+hwa2*a1x4*a1x4;        h45 += hwa1*(h2)+hwa2*a1x4*a1x5;
			h55 += hwa1*(h3)+hwa2*a1x5*a1x5;
		}
		h00 += hwb1*(ss1 + ss2) + hwb2*a2x0*a2x0; h01 += hwb1*(-ss1) + hwb2*a2x0*a2x1; h02 += hwb1*(-ss2) + hwb2*a2x0*a2x2; h03 += hwb2*a2x0*a2x3; h04 += hwb1*(-h4) + hwb2*a2x0*a2x4; h05 += hwb1*(h4)+hwb2*a2x0*a2x5;
		h11 += hwb1*(h1)+hwb2*a2x1*a2x1;        h12 += hwb1*(h2)+hwb2*a2x1*a2x2;    h13 += hwb1*(h4)+hwb2*a2x1*a2x3;    h14 += hwb2*a2x1*a2x4; h15 += hwb1*(-h4) + hwb2*a2x1*a2x5;
		h22 += hwb1*(h3)+hwb2*a2x2*a2x2;        h23 += hwb1*(-h4) + hwb2*a2x2*a2x3; h24 += hwb1*(h4)+hwb2*a2x2*a2x4;    h25 += hwb2*a2x2*a2x5;
		h33 += hwb1*(ss1 + ss2) + hwb2*a2x3*a2x3; h34 += hwb1*(-ss1) + hwb2*a2x3*a2x4; h35 += hwb1*(-ss2) + hwb2*a2x3*a2x5;
		h44 += hwb1*(h1)+hwb2*a2x4*a2x4;        h45 += hwb1*(h2)+hwb2*a2x4*a2x5;
		h55 += hwb1*(h3)+hwb2*a2x5*a2x5;

		u = aa*walpha; v = bb*wbeta;
		pardiso_b[f0] -= (u*a1x0 + v*a2x0);
		pardiso_b[f1] -= (u*a1x1 + v*a2x1);
		pardiso_b[f2] -= (u*a1x2 + v*a2x2);
		pardiso_b[f0 + V_N] -= (u*a1x3 + v*a2x3);
		pardiso_b[f1 + V_N] -= (u*a1x4 + v*a2x4);
		pardiso_b[f2 + V_N] -= (u*a1x5 + v*a2x5);

		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;
	}

	double len01, len0i, len1i, coef1, coef2, dis;
	int i;
	double s0, s1, e0, e1, p0, p1;
	for (int ii = 0; ii < AV_F_N; ++ii)
	{
		f0 = V_F0[AV_ID[ii]];
		f1 = V_F1[AV_ID[ii]];
		f2 = V_F2[AV_ID[ii]];
		x0 = position[f0];
		y0 = position[f0 + V_N];
		x1 = position[f1];
		y1 = position[f1 + V_N];
		x2 = position[f2];
		y2 = position[f2 + V_N];

		s0 = x0;	s1 = y0;
		e0 = x1;	e1 = y1;
		p0 = x2;	p1 = y2;
		len01 = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
		len0i = sqrt((p0 - s0)*(p0 - s0) + (p1 - s1)*(p1 - s1));
		len1i = sqrt((p0 - e0)*(p0 - e0) + (p1 - e1)*(p1 - e1));
		dis = len0i + len1i - len01;

		h_u = 2 * barrer_coef * threhold * (dis - threhold) / dis / dis / dis;
		uu = 2 * barrer_coef * threhold * (3 * threhold - 2 * dis) / dis / dis / dis / dis;
		a1x0 = (x0 - x2) / len0i - (x0 - x1) / len01;
		a1x1 = (x1 - x2) / len1i - (x1 - x0) / len01;
		a1x2 = (x2 - x0) / len0i + (x2 - x1) / len1i;
		a1x3 = (y0 - y2) / len0i - (y0 - y1) / len01;
		a1x4 = (y1 - y2) / len1i - (y1 - y0) / len01;
		a1x5 = (y2 - y0) / len0i + (y2 - y1) / len1i;

		coef1 = len01 * len01;
		coef2 = h_u / coef1 / len01;
		aa = (x0 - x1)*(x0 - x1);
		bb = (x1 - x0)*(y1 - y0);
		cc = (y0 - y1)*(y0 - y1);
		h00 = coef2 * (aa - coef1); h01 = coef2 * (coef1 - aa); h02 = 0.0; h03 = coef2 * bb; h04 = -coef2 * bb; h05 = 0.0;
		h11 = coef2 * (aa - coef1); h12 = 0.0; h13 = -coef2 * bb; h14 = coef2 * bb; h15 = 0.0;
		h22 = 0.0; h23 = 0.0; h24 = 0.0; h25 = 0.0;
		h33 = coef2 * (cc - coef1); h34 = coef2 * (coef1 - cc); h35 = 0.0;
		h44 = coef2 * (cc - coef1); h45 = 0.0;
		h55 = 0.0;

		h00 += uu*a1x0*a1x0; h01 += uu*a1x0*a1x1; h02 += uu*a1x0*a1x2; h03 += uu*a1x0*a1x3; h04 += uu*a1x0*a1x4; h05 += uu*a1x0*a1x5;
		h11 += uu*a1x1*a1x1; h12 += uu*a1x1*a1x2; h13 += uu*a1x1*a1x3; h14 += uu*a1x1*a1x4; h15 += uu*a1x1*a1x5;
		h22 += uu*a1x2*a1x2; h23 += uu*a1x2*a1x3; h24 += uu*a1x2*a1x4; h25 += uu*a1x2*a1x5;
		h33 += uu*a1x3*a1x3; h34 += uu*a1x3*a1x4; h35 += uu*a1x3*a1x5;
		h44 += uu*a1x4*a1x4; h45 += uu*a1x4*a1x5;
		h55 += uu*a1x5*a1x5;

		pardiso_b[f0] -= h_u * (a1x0);
		pardiso_b[f1] -= h_u * (a1x1);
		pardiso_b[f2] -= h_u * (a1x2);
		pardiso_b[f0 + V_N] -= h_u * (a1x3);
		pardiso_b[f1 + V_N] -= h_u * (a1x4);
		pardiso_b[f2 + V_N] -= h_u * (a1x5);

		i = F_N + AV_ID[ii];
		pardiso_a[id_h00[i]] += h00; pardiso_a[id_h01[i]] += h01; pardiso_a[id_h02[i]] += h02; pardiso_a[id_h03[i]] += h03; pardiso_a[id_h04[i]] += h04; pardiso_a[id_h05[i]] += h05;
		pardiso_a[id_h11[i]] += h11; pardiso_a[id_h12[i]] += h12; pardiso_a[id_h13[i]] += h13; pardiso_a[id_h14[i]] += h14; pardiso_a[id_h15[i]] += h15;
		pardiso_a[id_h22[i]] += h22; pardiso_a[id_h23[i]] += h23; pardiso_a[id_h24[i]] += h24; pardiso_a[id_h25[i]] += h25;
		pardiso_a[id_h33[i]] += h33; pardiso_a[id_h34[i]] += h34; pardiso_a[id_h35[i]] += h35;
		pardiso_a[id_h44[i]] += h44; pardiso_a[id_h45[i]] += h45;
		pardiso_a[id_h55[i]] += h55;
	}

	pardiso_a[0] += 1.0;
	pardiso_a[pardiso_ia[V_N]] += 1.0;
	pardiso->a = pardiso_a;
	pardiso->rhs = pardiso_b;
	long time_beg, time_end;
	time_beg = clock();
	pardiso->factorize();
	time_end = clock();
	double time_consumption = (time_end - time_beg) / 1000.0;
//	std::cout << "factorize" << time_consumption << std::endl;
	pardiso->pardiso_solver();

	vector<double> result_d = pardiso->result;
	VectorXd negative_grad(2 * V_N), d(2 * V_N);
	for (int i = 0; i < 2 * V_N; i++)
	{
		negative_grad(i) = pardiso_b[i];
		d(i) = result_d[i];
	}
	pardiso->free_numerical_factorization_memory();

	double temp_t;
	max_step(position_of_mesh, d, temp_t);
	double alpha = 0.95*temp_t;
	backtracking_line_search(position_of_mesh, d, negative_grad, alpha);
	position_of_mesh += alpha * d;
	Energysource();
}

void Parafun::max_step(const VectorXd &xx, const VectorXd &qq, double &step)
{
	double temp_t = numeric_limits<double>::infinity();
	int f0, f1, f2;
	double a, b, c, tt, tt1, tt2;
	double x0, x1, x2, y0, y1, y2, u0, u1, u2, v0, v1, v2;
	const double *x = xx.data();

	for (int i = 0; i < F_N; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = x[f0]; x1 = x[f1]; x2 = x[f2]; y0 = x[f0 + V_N]; y1 = x[f1 + V_N]; y2 = x[f2 + V_N];
		u0 = qq[f0]; u1 = qq[f1]; u2 = qq[f2]; v0 = qq[f0 + V_N]; v1 = qq[f1 + V_N]; v2 = qq[f2 + V_N];

		a = (u1 - u0) * (v2 - v0) - (v1 - v0) * (u2 - u0);
		b = (u1 - u0) * (y2 - y0) + (x1 - x0) * (v2 - v0) - (y1 - y0) * (u2 - u0) - (x2 - x0) * (v1 - v0);
		c = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
		tt = get_smallest_pos_quad_zero(a, b, c);

		if (temp_t > tt)
		{
			temp_t = tt;
		}
	}

	temp_t = 0.95*temp_t;
	const double *dir = qq.data();
	tmaxdetect(xx, qq, temp_t);
	double a2, b2, c2, a3, b3, c3;
	double root3[3];
	double d, q, r, r2, q3, A, B, tt3, tt4;
	double q1, q2, p1, p2, D, sqD, y, m2, m, n, alpha, beta, gamma, delta;
	double x0_up, x1_up, x2_up, y0_up, y1_up, y2_up;

	for (int i = 0; i < AV_F_N; ++i)
	{
		is_active[AV_ID[i]] = -1;

		f0 = V_F0[AV_ID[i]];
		f1 = V_F1[AV_ID[i]];
		f2 = V_F2[AV_ID[i]];
		x0 = x[f0]; x1 = x[f1]; x2 = x[f2];
		y0 = x[f0 + V_N]; y1 = x[f1 + V_N]; y2 = x[f2 + V_N];
		u0 = dir[f0]; u1 = dir[f1]; u2 = dir[f2];
		v0 = dir[f0 + V_N]; v1 = dir[f1 + V_N]; v2 = dir[f2 + V_N];

		//no triangle flip
		a = (u1 - u0) * (v2 - v0) - (v1 - v0) * (u2 - u0);
		b2 = (y1 - y0) * (u2 - u0) + (x2 - x0) * (v1 - v0);
		b = (u1 - u0) * (y2 - y0) + (x1 - x0) * (v2 - v0) - b2;
		c = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
		tt = numeric_limits<double>::infinity();
		if (b*b - 4 * a*c >= 0)
		{
			tt1 = 1 / (2 * a) * (-b + sqrt(b*b - 4 * a*c));
			tt2 = 1 / (2 * a) * (-b - sqrt(b*b - 4 * a*c));
			if (tt1 > 0 && tt2 > 0)
			{
				tt = min(tt1, tt2);
			}
			if (tt1 > 0 && tt2 < 0)
			{
				tt = tt1;
			}
			if (tt1 < 0 && tt2 > 0)
			{
				tt = tt2;
			}
		}

		if (temp_t > tt)
		{
			x0_up = x0 + tt * u0;	x1_up = x1 + tt * u1;	x2_up = x2 + tt * u2;
			y0_up = y0 + tt * v0;	y1_up = y1 + tt * v1;	y2_up = y2 + tt * v2;

			double dot1 = (x1_up - x0_up) * (x2_up - x0_up) + (y1_up - y0_up) * (y2_up - y0_up);
			double dot2 = (x0_up - x1_up) * (x2_up - x1_up) + (y0_up - y1_up) * (y2_up - y1_up);

			if (dot1 > 0 && dot2 > 0)
			{
				temp_t = tt;
			}
		}
	}

	step = temp_t;
	//std::cout << step << std::endl;
}
void Parafun::tmaxdetect(const VectorXd & x, const VectorXd &d, double& tmax)
{
	const double *pos = x.data();
	vector<double> x_update(2 * V_N);
	for (int i = 0; i < 2 * V_N; ++i)
	{
		x_update[i] = x[i] + tmax * d[i];
	}
	x_min = pos[0];
	x_max = pos[0];
	y_min = pos[V_N];
	y_max = pos[V_N];

	double x_min_up = x_update[0];
	double x_max_up = x_update[0];
	double y_min_up = x_update[V_N];
	double y_max_up = x_update[V_N];

	for (int i = 1; i < V_N; ++i)
	{
		if (pos[i] < x_min)
		{
			x_min = pos[i];
		}
		else if (pos[i] > x_max)
		{
			x_max = pos[i];
		}

		if (pos[i + V_N] < y_min)
		{
			y_min = pos[i + V_N];
		}
		else if (pos[i + V_N] > y_max)
		{
			y_max = pos[i + V_N];
		}

		if (x_update[i] < x_min_up)
		{
			x_min_up = x_update[i];
		}
		else if (x_update[i] > x_max_up)
		{
			x_max_up = x_update[i];
		}

		if (x_update[i + V_N] < y_min_up)
		{
			y_min_up = x_update[i + V_N];
		}
		else if (x_update[i + V_N] > y_max_up)
		{
			y_max_up = x_update[i + V_N];
		}
	}
	x_min = min(x_min, x_min_up);
	x_max = max(x_max, x_max_up);
	y_min = min(y_min, y_min_up);
	y_max = max(y_max, y_max_up);

	lengthgrid_x = (x_max - x_min) / (cellx_num - 1);
	x_max = x_min + cellx_num * lengthgrid_x;
	lengthgrid_y = (y_max - y_min) / (celly_num - 1);
	y_max = y_min + celly_num * lengthgrid_y;

	for (int j = 0; j < cell_points.size(); ++j)
	{
		cell_points[j].clear();
	}

	VectorXi boundary_vertex = d_.frame_ids;
	int id, id_x_min, id_x_max, id_y_min, id_y_max;
	double s0, s1, e0, e1, p0, p1, s0_, s1_, e0_, e1_;
	double l_x_min, l_x_max, l_y_min, l_y_max;

	AV_ID.clear();
	AV_F_N = 0;
	int id_start, id_end;
	Vector4d localx, localy;
	double a, len, b;

	for (int i = 0; i < BE_N; ++i)
	{
		id = boundary_vertex(i);
		s0 = pos[id];	s1 = pos[id + V_N];
		e0 = x_update[id]; e1 = x_update[id + V_N];

		l_x_min = min(s0, e0);
		l_x_max = max(s0, e0);
		l_y_min = min(s1, e1);
		l_y_max = max(s1, e1);

		id_x_min = floor((l_x_min - x_min) / lengthgrid_x);
		id_x_max = floor((l_x_max - x_min) / lengthgrid_x);
		id_y_min = floor((l_y_min - y_min) / lengthgrid_y);
		id_y_max = floor((l_y_max - y_min) / lengthgrid_y);
		if (id_y_max > celly_num - 1)
		{
			id_y_max = celly_num - 1;
		}
		if (id_x_max > cellx_num - 1)
		{
			id_x_max = cellx_num - 1;
		}

		for (int j_y = id_y_min; j_y < id_y_max + 1; ++j_y)
		{
			for (int j_x = id_x_min; j_x < id_x_max + 1; ++j_x)
			{
				int id_grid = j_x + j_y * cellx_num;
				cell_points[id_grid].push_back(id);
			}
		}
	}

	for (int i = 0; i < BE_N; ++i)
	{
		id_start = boundary_vertex(i);
		id_end = boundary_vertex((i + 1) % BE_N);
		s0 = pos[id_start];	s1 = pos[id_start + V_N];
		e0 = pos[id_end];	e1 = pos[id_end + V_N];
		s0_ = x_update[id_start];	s1_ = x_update[id_start + V_N];
		e0_ = x_update[id_end];	e1_ = x_update[id_end + V_N];

		localx[0] = s0;	localx[1] = e0;	localx[2] = s0_;	localx[3] = e0_;
		localy[0] = s1;	localy[1] = e1;	localy[2] = s1_;	localy[3] = e1_;
		l_x_min = localx.minCoeff();
		l_x_max = localx.maxCoeff();
		l_y_min = localy.minCoeff();
		l_y_max = localy.maxCoeff();

		len = sqrt((l_x_max - l_x_min)*(l_x_max - l_x_min) + (l_y_max - l_y_min)*(l_y_max - l_y_min));

		id_x_min = floor((l_x_min - x_min) / lengthgrid_x);
		id_x_max = floor((l_x_max - x_min) / lengthgrid_x);
		id_y_min = floor((l_y_min - y_min) / lengthgrid_y);
		id_y_max = floor((l_y_max - y_min) / lengthgrid_y);

		if (id_y_max > celly_num - 1)
		{
			id_y_max = celly_num - 1;
		}
		if (id_x_max > cellx_num - 1)
		{
			id_x_max = cellx_num - 1;
		}
		if (id_x_min < 0)
		{
			id_x_min = 0;
		}
		if (id_y_min < 0)
		{
			id_y_min = 0;
		}

		for (int j_y = id_y_min; j_y < id_y_max + 1; ++j_y)
		{
			for (int j_x = id_x_min; j_x < id_x_max + 1; ++j_x)
			{
				int id_grid = j_x + j_y * cellx_num;
				int size = cell_points[id_grid].size();
				for (int k = 0; k < size; ++k)
				{
					int id_mid = cell_points[id_grid][k];
					if (id_start == id_mid || id_end == id_mid)
					{
						continue;
					}

					if (BE_N - 1 == i)
					{
						if (is_active[i*(BE_N - 2) + boundary_vertexID[id_mid] - 1] < 0)
						{
							is_active[i*(BE_N - 2) + boundary_vertexID[id_mid] - 1] = 1;
							AV_ID.push_back(i*(BE_N - 2) + boundary_vertexID[id_mid] - 1);
						}
					}
					else if (boundary_vertexID[id_mid] < boundary_vertexID[id_start])
					{
						if (is_active[i*(BE_N - 2) + boundary_vertexID[id_mid]]  < 0)
						{
							is_active[i*(BE_N - 2) + boundary_vertexID[id_mid]] = 1;
							AV_ID.push_back(i*(BE_N - 2) + boundary_vertexID[id_mid]);
						}

					}
					else
					{
						if (is_active[i*(BE_N - 2) + boundary_vertexID[id_mid] - 2] < 0)
						{
							is_active[i*(BE_N - 2) + boundary_vertexID[id_mid] - 2] = 1;
							AV_ID.push_back(i*(BE_N - 2) + boundary_vertexID[id_mid] - 2);
						}
					}
				}
			}
		}
	}

	AV_F_N = AV_ID.size();
	//	std::cout << AV_F_N << std::endl;
}
double Parafun::get_smallest_pos_quad_zero(double a, double b, double c)
{
	using namespace std;
	double t1, t2;
	if (std::abs(a) < 1.0e-10)
	{
		a *= 1e6;
		b *= 1e6;
		c *= 1e6;
	}
	if (std::abs(a) > 1.0e-10)
	{
		double delta_in = pow(b, 2) - 4 * a * c;
		if (delta_in <= 0)
		{
			return INFINITY;
		}

		double delta = sqrt(delta_in); // delta >= 0
		if (b >= 0) // avoid subtracting two similar numbers
		{
			double bd = -b - delta;
			t1 = 2 * c / bd;
			t2 = bd / (2 * a);
		}
		else
		{
			double bd = -b + delta;
			t1 = bd / (2 * a);
			t2 = (2 * c) / bd;
		}

		assert(std::isfinite(t1));
		assert(std::isfinite(t2));

		if (a < 0) std::swap(t1, t2); // make t1 > t2
									  // return the smaller positive root if it exists, otherwise return infinity
		if (t1 > 0)
		{
			return t2 > 0 ? t2 : t1;
		}
		else
		{
			return INFINITY;
		}
	}
	else
	{
		if (b == 0) return INFINITY; // just to avoid divide-by-zero
		t1 = -c / b;
		return t1 > 0 ? t1 : INFINITY;
	}
}

void Parafun::backtracking_line_search(const VectorXd &x, const VectorXd &d, const VectorXd &negetive_grad, double &alpha, bool is_interp)
{
	double h = 0.5;
	double tt = -(negetive_grad.transpose()*d)(0, 0);
	double c = 0.2;
	double ex;
	Energy(x, ex, is_interp, false);
	ex = ex + barrer_coef * energy_barrier;
	double e;
	VectorXd x_new = x + alpha * d;
	Energy(x_new, e, is_interp);
	while (e > ex + alpha * c * tt)
	{
		alpha = h*alpha;
		x_new = x + alpha * d;
		Energy(x_new, e, is_interp);
	}
}
void Parafun::Energy(const VectorXd &position, double &energyupdate, bool is_interp, bool is_whole)
{
	double energy = 0;

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_d;
	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;
	const double *pos = position.data();

	double* tmp_p00;
	double* tmp_p01;
	double* tmp_p10;
	double* tmp_p11;

	if (is_interp)
	{
		tmp_p00 = update_p00.data();
		tmp_p01 = update_p01.data();
		tmp_p10 = update_p10.data();
		tmp_p11 = update_p11.data();
	}
	else
	{
		tmp_p00 = source_p00.data();
		tmp_p01 = source_p01.data();
		tmp_p10 = source_p10.data();
		tmp_p11 = source_p11.data();
	}

	for (int i = 0; i < F_N; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + total_num];

		x1 = pos[f1];
		y1 = pos[f1 + total_num];

		x2 = pos[f2];
		y2 = pos[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = tmp_p00[i]; p01 = tmp_p01[i]; p10 = tmp_p10[i]; p11 = tmp_p11[i];

		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;


		det = j00*j11 - j01*j10;
		E_d = (1 + 1 / (det*det)) * (j00*j00 + j01*j01 + j10*j10 + j11*j11);

		energy += area[i] * E_d;
	}

	energyupdate = energy;

	if (is_whole)
	{
		fungrid(position);
		double dis, E_b, energy2 = 0;
		for (int i = 0; i < AV_F_N; ++i)
		{
			f0 = V_F0[AV_ID[i]];
			f1 = V_F1[AV_ID[i]];
			f2 = V_F2[AV_ID[i]];

			x0 = pos[f0];	y0 = pos[f0 + V_N];

			x1 = pos[f1];	y1 = pos[f1 + V_N];

			x2 = pos[f2];	y2 = pos[f2 + V_N];

			dis = distance(x0, y0, x1, y1, x2, y2);
			if (dis < 0)
			{
				std::cout << "distance is zero" << std::endl;
			}
			E_b = (1 - threhold / dis) * (1 - threhold / dis);
			energy2 += E_b;
		}
		energyupdate = energy + barrer_coef * energy2;
		energy_barrier = energy2;
		//std::cout << AV_F_N << std::endl;
	}
}

void Parafun::Energysource()
{
	double end_e_area = 0;

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_1, E_2;

	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	const double *pos = position_of_mesh.data();
	for (int i = 0; i < d_.m_T.rows(); ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + total_num];

		x1 = pos[f1];
		y1 = pos[f1 + total_num];

		x2 = pos[f2];
		y2 = pos[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;

		det = j00*j11 - j01*j10;

		E_1 = (j00*j00 + j01*j01 + j10*j10 + j11*j11);
		E_2 = 1.0 / (det*det)* E_1;

		end_e_area += ((E_1 + E_2)*area[i]);
	}

	energy_mesh = end_e_area;

	end_e_area = 0;
	for (int i = d_.m_T.rows(); i < F_N; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + total_num];

		x1 = pos[f1];
		y1 = pos[f1 + total_num];

		x2 = pos[f2];
		y2 = pos[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00*q00 + p10*q01; j01 = p01*q00 + p11*q01; j10 = p00*q10 + p10*q11; j11 = p01*q10 + p11*q11;

		det = j00*j11 - j01*j10;

		E_1 = (j00*j00 + j01*j01 + j10*j10 + j11*j11);
		E_2 = 1.0 / (det*det)* E_1;

		end_e_area += (E_1 + E_2);
	}

	energy_shell = end_e_area;
}
double Parafun::compute_energy(const Eigen::MatrixXd & x, bool whole)
{
	double end_e_one_temp = 0, end_e_area = 0;

	int f0, f1, f2;
	double x0, y0, x1, y1, x2, y2;
	double det, E_1, E_2;

	double j00, j01, j10, j11;
	double p00, p01, p10, p11;
	double q00, q01, q10, q11;

	const double *pos = x.data();
	int src_t_num = d_.m_T.rows();

	for (int i = 0; i < src_t_num; ++i)
	{
		f0 = F0[i];
		f1 = F1[i];
		f2 = F2[i];

		x0 = pos[f0];
		y0 = pos[f0 + total_num];

		x1 = pos[f1];
		y1 = pos[f1 + total_num];

		x2 = pos[f2];
		y2 = pos[f2 + total_num];

		q00 = x1 - x0; q01 = x2 - x0;
		q10 = y1 - y0; q11 = y2 - y0;

		p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

		j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

		det = j00 * j11 - j01 * j10;


		E_1 = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
		E_2 = 1.0 / (det*det)* E_1;

		end_e_one_temp += E_1;
		end_e_one_temp += E_2;
		end_e_area += ((E_1 + E_2)*area_src[i]);

	}

	if (whole)
	{
		for (int i = src_t_num; i < F_N; ++i)
		{
			f0 = F0[i];
			f1 = F1[i];
			f2 = F2[i];

			x0 = pos[f0];
			y0 = pos[f0 + total_num];

			x1 = pos[f1];
			y1 = pos[f1 + total_num];

			x2 = pos[f2];
			y2 = pos[f2 + total_num];

			q00 = x1 - x0; q01 = x2 - x0;
			q10 = y1 - y0; q11 = y2 - y0;

			p00 = source_p00[i]; p01 = source_p01[i]; p10 = source_p10[i]; p11 = source_p11[i];

			j00 = p00 * q00 + p10 * q01; j01 = p01 * q00 + p11 * q01; j10 = p00 * q10 + p10 * q11; j11 = p01 * q10 + p11 * q11;

			det = j00 * j11 - j01 * j10;

			E_1 = (j00*j00 + j01 * j01 + j10 * j10 + j11 * j11);
			E_2 = 1.0 / (det*det)* E_1;

			end_e_one_temp += E_1;
			end_e_one_temp += E_2;
			end_e_area += ((E_1 + E_2)*area[i]);
		}

		end_e_area += barrer_coef * energy_barrier;
	}
	//cout << "compute energy with scaf "<<whole<<" "<< end_e_area << endl;
	//cout << "compute energy with scaf/area " << whole << " " << end_e_area/d_.mesh_measure << endl;

	return end_e_area;
}

void Parafun::local_coordinate_inverse(int i, double &p00, double &p01, double &p10, double &p11)
{
	int f0 = F0[i];
	int f1 = F1[i];
	int f2 = F2[i];

	Vector3d x_(d_.m_V(f1, 0) - d_.m_V(f0, 0), d_.m_V(f1, 1) - d_.m_V(f0, 1), d_.m_V(f1, 2) - d_.m_V(f0, 2));
	double x1_0 = x_.norm();
	x_ /= x1_0;
	Vector3d l_(d_.m_V(f2, 0) - d_.m_V(f0, 0), d_.m_V(f2, 1) - d_.m_V(f0, 1), d_.m_V(f2, 2) - d_.m_V(f0, 2));

	Vector3d n_ = x_.cross(l_);
	n_.normalize();
	Vector3d y_ = n_.cross(x_);
	double x2_0 = l_.dot(x_);
	double y2_0 = l_.dot(y_);

	p00 = 1 / x1_0;
	p01 = -x2_0 / (x1_0*y2_0);
	p10 = 0;
	p11 = 1 / y2_0;
}
void Parafun::local_coordinate_inverse_scaf(int i, double & p00, double & p01, double & p10, double & p11)
{
	int f0 = F0[i];
	int f1 = F1[i];
	int f2 = F2[i];

	Vector2d x_(d_.w_uv(f1, 0) - d_.w_uv(f0, 0), d_.w_uv(f1, 1) - d_.w_uv(f0, 1));
	Vector2d l_(d_.w_uv(f2, 0) - d_.w_uv(f0, 0), d_.w_uv(f2, 1) - d_.w_uv(f0, 1));

	double area_tri = abs(x_(0)*l_(1) - x_(1)*l_(0));
	double x1_0, x2_0, y2_0;
	if (area_tri > area_threshold)
	{
		x1_0 = x_.norm();
		x_ /= x1_0;
		Vector2d y_(-x_(1), x_(0));
		x2_0 = l_.dot(x_);
		y2_0 = l_.dot(y_);
	}
	else
	{
		//cout << "area too small!!!!!!!!!!!!! " << endl;
		double h = sqrt((2 * area_threshold) / sqrt(3.0));
		x1_0 = h;
		x2_0 = h / 2.0;
		y2_0 = sqrt(3.0)*h / 2.0;
	}
	p00 = 1 / x1_0;
	p01 = -x2_0 / (x1_0*y2_0);
	p10 = 0;
	p11 = 1 / y2_0;
}

double Parafun::newton_equation(const double & a, const double & b, const double & K)
{
	double tt = 1;
	double E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K;
	while (abs(E_d) > 1e-5)
	{
		tt = tt - 1 / (2 * log(a)*pow(a, 2 * tt) + 2 * log(b)* pow(b, 2 * tt) + 2 * log(1 / a)* pow(1 / a, 2 * tt) + 2 * log(1 / b)*pow(1 / b, 2 * tt))*(pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K);
		E_d = pow(a, 2 * tt) + pow(b, 2 * tt) + pow(1 / a, 2 * tt) + pow(1 / b, 2 * tt) - K;
	}
	return tt;
}
void Parafun::adjust_shell_weight(double new_weight)
{
	d_.shell_factor = new_weight;
	d_.update_shell();
	init_area();
}

double Parafun::distance(double s0, double s1, double e0, double e1, double p0, double p1)
{
	double s0_ = s0, s1_ = s1;
	double e0_ = e0, e1_ = e1;

	double len01 = sqrt((s0 - e0)*(s0 - e0) + (s1 - e1)*(s1 - e1));
	double len0i = sqrt((p0 - s0_)*(p0 - s0_) + (p1 - s1_)*(p1 - s1_));
	double len1i = sqrt((p0 - e0_)*(p0 - e0_) + (p1 - e1_)*(p1 - e1_));

	return len0i + len1i - len01;
}

void Parafun::fungrid(const VectorXd &x)
{
	const double *pos = x.data();

	x_min = pos[0];
	x_max = pos[0];
	y_min = pos[V_N];
	y_max = pos[V_N];
	for (int i = 1; i < V_N; ++i)
	{
		if (pos[i] < x_min)
		{
			x_min = pos[i];
		}
		else if (pos[i] > x_max)
		{
			x_max = pos[i];
		}

		if (pos[i + V_N] < y_min)
		{
			y_min = pos[i + V_N];
		}
		else if (pos[i + V_N] > y_max)
		{
			y_max = pos[i + V_N];
		}
	}

	lengthgrid_x = (x_max - x_min) / (cellx_num - 1);
	x_max = x_min + cellx_num * lengthgrid_x;
	lengthgrid_y = (y_max - y_min) / (celly_num - 1);
	y_max = y_min + celly_num * lengthgrid_y;

	for (int j = 0; j < cell_points.size(); ++j)
	{
		cell_points[j].clear();
	}
	VectorXi boundary_vertex = d_.frame_ids;
	int bound_num = BE_N;
	double l_x_min, l_x_max, l_y_min, l_y_max, b, len;
	int id_x_min, id_x_max, id_y_min, id_y_max;
	int id_start, id_end;
	double s0, s1, e0, e1, p0, p1, s0_, s1_, e0_, e1_;

	AV_ID.clear();
	AV_F_N = 0;
	double dis;

	for (int i = 0; i < bound_num; ++i)
	{
		int id = boundary_vertex(i);
		int x_i = std::floor((pos[id] - x_min) / lengthgrid_x);
		int y_i = std::floor((pos[id + V_N] - y_min) / lengthgrid_y);
		if (y_i > celly_num - 1)
		{
			y_i = celly_num - 1;
		}
		if (x_i > cellx_num - 1)
		{
			x_i = cellx_num - 1;
		}
		cell_points[y_i * cellx_num + x_i].push_back(id);
	}

	for (int i = 0; i < bound_num; ++i)
	{
		id_start = boundary_vertex(i);
		id_end = boundary_vertex((i + 1) % BE_N);
		s0 = pos[id_start];	s1 = pos[id_start + V_N];
		e0 = pos[id_end];	e1 = pos[id_end + V_N];
		len = sqrt((s0 - e0)*(s0 - e0) + (s1 - e1)*(s1 - e1));
		b = sqrt(threhold*len / 2 + threhold*threhold / 4);
		s0_ = s0;
		s1_ = s1;
		e0_ = e0;
		e1_ = e1;
		l_x_min = min(s0_, e0_);
		l_x_max = max(s0_, e0_);
		l_y_min = min(s1_, e1_);
		l_y_max = max(s1_, e1_);

		//local box
		l_x_min = max(l_x_min - b, x_min);
		l_x_max = min(l_x_max + b, x_max);
		l_y_min = max(l_y_min - b, y_min);
		l_y_max = min(l_y_max + b, y_max);

		id_x_min = floor((l_x_min - x_min) / lengthgrid_x);
		id_x_max = floor((l_x_max - x_min) / lengthgrid_x);
		id_y_min = floor((l_y_min - y_min) / lengthgrid_y);
		id_y_max = floor((l_y_max - y_min) / lengthgrid_y);

		if (id_y_max > celly_num - 1)
		{
			id_y_max = celly_num - 1;
		}
		if (id_x_max > cellx_num - 1)
		{
			id_x_max = cellx_num - 1;
		}

		for (int j_y = id_y_min; j_y < id_y_max + 1; ++j_y)
		{
			for (int j_x = id_x_min; j_x < id_x_max + 1; ++j_x)
			{
				int id_cell = j_x + j_y * cellx_num;
				int size = cell_points[id_cell].size();
				for (int k = 0; k < size; ++k)
				{
					int id_p = cell_points[id_cell][k];

					if (id_start == id_p || id_end == id_p)
					{
						continue;
					}

					p0 = pos[id_p];	p1 = pos[id_p + V_N];
					dis = sqrt((p0 - s0_)*(p0 - s0_) + (p1 - s1_)*(p1 - s1_)) + sqrt((p0 - e0_)*(p0 - e0_) + (p1 - e1_)*(p1 - e1_)) - len;
					if (dis > threhold)
					{
						continue;
					}
					if (dis < 0)
					{
						std::cout << "distance error" << std::endl;
					}

					if (BE_N - 1 == i)
					{
						AV_ID.push_back(i*(BE_N - 2) + boundary_vertexID[id_p] - 1);
					}
					else if (id_p < id_start)
					{
						AV_ID.push_back(i*(BE_N - 2) + boundary_vertexID[id_p]);
					}
					else
					{
						AV_ID.push_back(i*(BE_N - 2) + boundary_vertexID[id_p] - 2);
					}
				}
			}
		}
	}

	AV_F_N = AV_ID.size();
	//std::cout << AV_F_N << std::endl;
}

Parafun::~Parafun()
{
	if (pardiso != NULL)
	{
		delete pardiso;
		pardiso = NULL;
	}
}
