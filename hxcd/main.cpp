#include <iostream>
#include <vector>
#include <array>
#include "Transport.h"
#include "UserSetting.h"
#include <sstream>
#include <iomanip>
#include "time.h"
using namespace std;

static char help[] = "Solve Motor-Assisted Transport Equation\n";

int main(int argc, char **argv)
{
	if (argc == 3)
	{
		stringstream ss, stmp;
		string path_in;
		ss << argv[1];
		ss >> path_in;
		stmp << argv[2];
		int n_process = atoi(argv[2]);

		int rank, nProcs;
		PetscErrorCode ierr;
		/// start up petsc
		ierr = PetscInitialize(&argc, &argv, (char*)0, help); if (ierr) return ierr;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
		MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);

		int n_bzmesh;
		vector<double> var;
		vector<Vertex2D> cpts;
		vector<array<double, 2>> velocity_node;
		vector<vector<int>> ele_process;
		vector<Element2D> tmesh;
		vector<double>  N0_ini, Nplus_ini, Nminus_ini;
		ele_process.resize(nProcs);
		

		/// Set simulation parameters and mesh
		string fn_mesh(path_in + "controlmesh.vtk");
		string fn_bz(path_in + "bzmeshinfo.txt.epart." + stmp.str());
		string fn_velocity(path_in + "velocityfield.txt");		
		string fn_parameter(path_in+"simulation_parameter.txt");
		string path_out(path_in);
		
		UserSetting *user = new UserSetting;
		user->ReadMesh(fn_mesh, cpts, tmesh);
		user->ReadVelocityField(fn_velocity, cpts.size(), velocity_node);
		user->AssignProcessor(fn_bz, n_bzmesh, ele_process);
		user->SetVariables(fn_parameter, var, path_out);
		user->SetInitialCondition(cpts, var, N0_ini, Nplus_ini, Nminus_ini);

		///Heat Convection-diffusion 2D steady Problem
		Transport* transport = new Transport;
		transport->InitializeProblem(n_bzmesh, velocity_node, N0_ini, Nplus_ini, Nminus_ini, var); //CA = CA0
		transport->AssignProcessor(ele_process);
		transport->Run(cpts, velocity_node, tmesh, path_in, path_out);

		PetscPrintf(PETSC_COMM_WORLD, "Done!\n");
		ierr = PetscFinalize(); CHKERRQ(ierr);
	}
	else if (argc > 3) {
		cout << "Too many arguments.\n";
	}
	else {
		cout << "Two argument expected.\n";
	}
	return 0;
}
