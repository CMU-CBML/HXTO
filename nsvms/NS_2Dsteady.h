#ifndef NS_2DSTEADY_H
#define NS_2DSTEADY_H

#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include "BasicDataStructure.h"
#include "time.h"

using namespace std;

const int degree = 3;
const int dim = 2;
const int bzpt_num = 16;

class NS_2Dsteady
{
private:
	vector<double> Gpt;
	vector<double> wght;
	const double PI = 4 * atan(1.0);

	PetscErrorCode ierr;
	MPI_Comm comm;
	int mpiErr;
	int comRank;
	int comSize;
	int nProcess;

	int rstart, rend;
	int n_bzmesh;
	vector<int> ele_process;
	vector<Element2D> bzmesh_process;

	KSP ksp;
	PC pc;
	Mat GK;
	Vec GR;
	Vec temp_solution;
	
	double dt;
	double velocity_max;
	double nu;
	double rou;
	double alphaM;
	double alphaF;
	double Gama;
	double fx, fy; //fz;//elemental body force

	vector<double> Vel;
	vector<double> Pre;
	vector<double> gh; //prescribed boundary velocities
	vector<double> par;//Dn0, v_plus, v_minus, k+, k-,k'+,k'-
	

public:
	NS_2Dsteady();
private:
	void BodyForce(double x, double y, double &Fx, double &Fy);
	void ReadBezierElementProcess(string fn);

	/*Analysis*/
	void GaussInfo(int ng);	
	void BasisFunction(double u, double v, int nen, const vector<array<double, 2>>& pt, const vector<array<double, 16>> &cmat, vector<double> &Nx, vector<array<double, 2>> &dNdx, vector<array<array<double, 2>, 2>> &dN2dx2, double dudx[2][2], double& detJ);
	void PointFormValue(vector<double>& Nx, const vector<array<double, 3>>& U, double Value[3]);
	void PointFormGrad(vector<array<double, 2>>& dNdx, const vector<array<double, 3>>& U, double Value[3][2]);
	void PointFormHess(vector<array<array<double, 2>, 2>>& d2Ndx2, const vector<array<double, 3>> &U, double Value[3][2][2]);
	void Tau(double J[2][2], double u[3], double &tauM, double &tauC);
	void FineScale(double tauM, double tauC, double u[2], double u_x[2], double u_y[2], double u_xx[2], double u_yy[2], double p, double p_x, double p_y, double u_s[2], double &p_s);
	void Residual(vector<double>& Nx, vector<array<double, 2>>& dNdx, vector<array<array<double, 2>, 2>>& dN2dx2, double dudx[2][2], const double detJ, const vector<array<double, 3>> &U, vector<array<double, 3>> Re);
	void Tangent(vector<double> &Nx, vector<array<double, 2>>& dNdx, double dudx[2][2], const double detJ, const vector<array<double, 3>>& U, vector<array<vector<array<double, 3>>, 3>>& Ke);
	void BuildLinearSystemProcess(const vector<Vertex2D>& cpts, const vector<array<double, 2>>& velocity_bc, const vector<double> velocity_node, const vector<double> pressure_node);
	void ApplyBoundaryCondition(const double bc_value, int pt_num, int variable_num, vector<array<vector<array<double, 3>>, 3>>& Ke, vector<array<double, 3>> &Re);
	void MatrixAssembly(vector<array<vector<array<double, 3>>, 3>>& Ke, const vector<int>& IEN, Mat& GK);
	void ResidualAssembly(vector<array<double, 3>> &Re, const vector<int>& IEN, Vec& GR);
  void Pressure_Drop(const vector<Vertex2D>& cpts, const vector<double> pressure_node, double& pre_drop);

	/*Postprocessing*/
	void ResultCal_Bezier(double u, double v, const Element2D& bzel, double pt[3], double result[3], double dudx[2], double& detJ);
	void VisualizeVTK_ControlMesh(const vector<Vertex2D>& pts, const vector<Element2D>& mesh, int step, string fn);
	void VisualizeVTK_PhysicalDomain(int step, string fn);
	void WriteVTK(const vector<array<double, 3>> pts, const vector<double> sdisp, const vector<array<int, 4>> sele, int step, string fn);
public:
	/*Preprocessing*/
	void InitializeProblem(const int ndof, const int n_bz, const vector<double>& Vel0, const vector<double>& Pre0, const vector<double>& var);
	void AssignProcessor(vector<vector<int>> &ele_proc);
	void Run(const vector<Vertex2D>& cpts, const vector<Element2D>& tmesh, const vector<array<double, 2>>& velocity_bc, string fn);
};

#endif
