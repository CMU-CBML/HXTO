#ifndef TRANSPORT_H
#define TRANSPORT_H

#include <vector>
#include <array>
#include "BasicDataStructure.h"
#include "time.h"

using namespace std;

class Transport
{
private:

	PetscErrorCode ierr;
	MPI_Comm comm;
	int mpiErr;
	int comRank;
	int comSize;
	int nProcess;
	int n_bzmesh;
	int rstart, rend;
	vector<int> ele_process;
	vector<Element2D> bzmesh_process;

	KSP ksp;
	PC pc;
	Mat GK;
	Vec GR;
	Vec temp_solution;

	int judge;// judge if the matrix have been solved
	vector<double> Gpt;
	vector<double> wght;
	vector<double> N_0;
	vector<double> N_plus;
	vector<double> N_minus;
	vector<array<double, 2>> vplus, vminus;
	vector<double> par;//Dn0, v_plus, v_minus, k+, k-,k'+,k'-, dt, nstep, n0_bc, n+_bc, n-_bc
	double dt;
	int nstep;
public:
	Transport();
private:
	void ReadBezierElementProcess(string fn);

	void GaussInfo(int ng);
	void SUPGcoefficient(double s, double v[2], double dudx[2][2], vector<array<double, 2>>& dNdx, double &tau_supg);
	void BasisFunction(double u, double v, const vector<array<double, 2>>& pt, double Nx[bzpt_num], double dNdx[bzpt_num][dim], double dudx[2][2], double& detJ);
	void BasisFunction(double u, double v, const vector<array<double, 2>>& pt, const vector<array<double, 16>> &cmat, vector<double> &Nx, vector<array<double, 2>> &dNdx, double dudx[2][2], double& detJ);
	void WeightingFunction(const double velocity[2], const double& s, const double& tau, const vector<double> &Nx, const vector<array<double, 2>> &dNdx, vector<double> &Wx);
	void ElementValue(const vector<double> &Nx, const vector<double> value_node, double &value);
	void ElementVelocity(const vector<double> &Nx, const vector<array<double, 2>>& v_node, double v_tmp[2]);
	void Tangent(const int nen, vector<double>& Nx, vector<double>& Npx, vector<array<double, 2>>& dNdx, double vp[2], double detJ, vector<vector<double>>& EMatrixSolve);
	void Residual(const int nen, const double CA, const double Nplus, const vector<double> &Nx, const vector<double> &Npx, const double detJ, vector<double> &EVectorSolve);
	void ApplyBoundaryCondition(const double bc_value, int pt_num, int variable_num, vector<vector<double>>& EMatrixSolve, vector<double>& EVectorSolve);
	void MatrixAssembly(vector<vector<double>>& EMatrixSolve, const vector<int>& IEN, Mat& GK);
	void ResidualAssembly(vector<double>& EVectorSolve, const vector<int>& IEN, Vec& GR);
	void BuildLinearSystemProcess(const vector<Element2D>& tmesh, const vector<Vertex2D> &cpts, const vector<array<double, 2>> velocity_node, const double Vplus, const double Vminus);

	void HeatExchange(const vector<Element2D>& tmesh, const vector<Vertex2D> &cpts, const vector<array<double, 2>> velocity_node, const double Vplus, vector<double> &heat_temp, double& heat_exchange);
	void ConcentrationCal_Coupling_Bezier(double u, double v, const Element2D& bzel, double pt[2], double& disp, double dudx[2], double& detJ);
	void VisualizeVTK_ControlMesh(const vector<Vertex2D>& pts, const vector<Element2D>& mesh, int step, string fn);
	void VisualizeVTK_PhysicalDomain(int step, string fn);
	void WriteVTK(const vector<array<double, 3>> pts, const vector<double> sdisp, const vector<array<int, 4>> sele, int step, string fn);
public:	
	void InitializeProblem(const int n_bz, vector<array<double, 2>> &velocity_node, const vector<double>& N0_ini, const vector<double>& Nplus_ini, const vector<double>& Nminus_ini, const vector<double>& var);
	void AssignProcessor(vector<vector<int>> &ele_proc);
	void Run(const vector<Vertex2D>& cpts, const vector<array<double, 2>> velocity_node, const vector<Element2D> &tmesh, string path_in, string path_out);

};

#endif