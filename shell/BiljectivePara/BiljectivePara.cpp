#include "BiljectivePara.h"

BiljectivePara::BiljectivePara(string filename)
{
	OpenMesh::IO::read_mesh(mesh, filename);
	string file_str = filename;
	modelname = file_str.substr(file_str.find_last_of('/') + 1);
	modelname.replace(modelname.end() - 4, modelname.end(), "");
	is_initialization = false;
	weight1 = true;
	weight2 = true;

	//update mesh center
	typedef Mesh::Point Point;
	Mesh::VertexIter vIt = mesh.vertices_begin();
	Mesh::VertexIter vEnd = mesh.vertices_end();
	OpenMesh::Vec3d bbMin = OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh.point(vIt));
	OpenMesh::Vec3d bbMax = OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh.point(vIt));

	size_t count = 0;
	for (; vIt != vEnd; ++vIt, ++count)
	{
		bbMin.minimize(OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh.point(vIt)));
		bbMax.maximize(OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh.point(vIt)));
	}

	OpenMesh::Vec3d c = (bbMin + bbMax)*0.5;
	for (vIt = mesh.vertices_begin(); vIt != vEnd; ++vIt, ++count)
	{
		OpenMesh::Vec3d p = mesh.point(vIt) - c;
		mesh.set_point(vIt, p);
	}
}

BiljectivePara::~BiljectivePara()
{
	
}

void BiljectivePara::parameterization()
{
	if (!is_initialization)
	{
		load();
	}

	long time_beg, time_end;
	time_beg = clock();
	int iteration_count = 0;
	double conv_rate_mesh = 1.0;
	double conv_rate_all = 1.0;
	bool is_first = true;
	bool is_second = true;

	double last_mesh_energy = parafun_solver->compute_energy(shell_data.w_uv, false) / shell_data.mesh_measure;
	parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure / (shell_data.sf_num) / 1000.0);
	
	for (int i = 0; i < MAX_ITER_NUM; i++) 
	{
		iteration_count++;
		
		bool is_ip_convrate = true;
		if (conv_rate_mesh < 0.01)
			is_ip_convrate = false;
		bool is_slim_convrate = false;
		if (conv_rate_mesh > 0.1)
			is_slim_convrate = true;

		shell_data.energy = parafun_solver->BPE(is_ip_convrate, is_slim_convrate);
		double current_mesh_energy = parafun_solver->energy_mesh / shell_data.mesh_measure;
		double mesh_energy_decrease = last_mesh_energy - current_mesh_energy;
		conv_rate_mesh = abs(mesh_energy_decrease) / last_mesh_energy;
		last_mesh_energy = current_mesh_energy;
		adjust_weight(conv_rate_mesh, last_mesh_energy);

		if (conv_rate_mesh < convgence_con_rate)
		{
			break;
		}
	}

	time_end = clock();
	double time_consumption = (time_end - time_beg) / 1000.0;
	
	cout << "COMP ====== time_consumption: " << time_consumption << ";	mesh energy: " << last_mesh_energy/4 << ";	sum_iter: " << iteration_count << endl;
	cout << modelname << " bijective parameterization finish !" << endl;

	write_obj();
}

double BiljectivePara::adjust_weight(double conv_mesh, double last_mesh_energy)
{
	if (conv_mesh > 1e-3 && weight1)
	{
		parafun_solver->adjust_shell_weight((last_mesh_energy)*shell_data.mesh_measure / (shell_data.sf_num) / 1000.0);
	}
	else if (conv_mesh > 1e-5 && weight2)
	{
		parafun_solver->adjust_shell_weight(4.0*shell_data.mesh_measure / (shell_data.sf_num) / 100000.0);
		weight1 = false;
	}
	else
	{
		parafun_solver->adjust_shell_weight(4.0*shell_data.mesh_measure / shell_data.sf_num / 100000000.0);
		weight2 = false;
	}

	double last_all_energy = parafun_solver->energy_mesh + parafun_solver->area.back() * parafun_solver->energy_shell + parafun_solver->barrer_coef * parafun_solver->energy_barrier;
	return last_all_energy;
}

void BiljectivePara::load()
{
	is_initialization = true;
	cout << modelname << " bijective parameterization begin ..." << endl;

	int v_num = mesh.n_vertices();
	int f_num = mesh.n_faces();
	Eigen::MatrixXd V(v_num, 3);
	Eigen::MatrixXi F(f_num, 3);

	for (int i = 0; i < v_num; ++i)
	{
		auto vertex = mesh.vertex_handle(i);
		Mesh::Point pos = mesh.point(vertex);
		for (int j = 0; j < 3; ++j)
		{
			V(i, j) = pos[j];
		}
	}

	for (int i = 0; i < f_num; ++i)
	{
		auto face = mesh.face_handle(i);
		int dd = 0;
		for (Mesh::FaceVertexIter it2 = mesh.fv_begin(face); it2 != mesh.fv_end(face); ++it2)
		{
			auto vertex_ = it2.handle();
			switch (dd)
			{
			case 0: F(i, 0) = vertex_.idx(); break;
			case 1: F(i, 1) = vertex_.idx(); break;
			case 2: F(i, 2) = vertex_.idx(); break;
			default:
				break;
			}
			dd++;
		}
	}
	
	shell_data = ShellData();
	shell_data.add_new_patch(V, F, Eigen::RowVector2d(0, 0));
	parafun_solver.reset(new Parafun(shell_data));
	
	for (int i = 0; i < v_num; ++i)
	{
		auto vertex = mesh.vertex_handle(i);
		auto p = shell_data.w_uv.row(i);
		Mesh::Point pos(p(0), p(1), 0.0);
		mesh.set_point(vertex, pos);
	}

	parafun_solver->after_mesh_improve();
}

void BiljectivePara::write_obj()
{
	Eigen::MatrixXd& V_in = shell_data.w_uv;
	int v_num = mesh.n_vertices();
	for (int i = 0; i < v_num; ++i)
	{
		auto vertex = mesh.vertex_handle(i);
		Mesh::Point pos(V_in(i, 0), V_in(i, 1), 0.0);
		mesh.set_point(vertex, pos);
	}

	Eigen::MatrixXi& F_ref = shell_data.m_T;
	string outstr = modelname + "_bijective_para_result.obj";

	ofstream of_obj(outstr, ios::trunc);

	if (V_in.cols() == 3)
	{
		for (size_t vid = 0; vid < V_in.rows(); vid++)
		{
			of_obj << "v " << V_in(vid, 0) << " " << V_in(vid, 1) << " " << V_in(vid, 2) << endl;
		}
	}
	else if (V_in.cols() == 2)
	{
		for (size_t vid = 0; vid < V_in.rows(); vid++)
		{
			of_obj << "v " << V_in(vid, 0) << " " << V_in(vid, 1) << " " << 0.0 << endl;
		}
	}

	for (size_t fi = 0; fi < F_ref.rows(); fi++)
	{
		of_obj << "f " << F_ref(fi, 0) + 1 << " " << F_ref(fi, 1) + 1 << " " << F_ref(fi, 2) + 1 << endl;
	}
	of_obj.close();
}

void BiljectivePara::shelltri(MatrixXi & tri, MatrixXd& pos)
{
	tri = shell_data.s_T;
	pos = shell_data.w_uv;
}
