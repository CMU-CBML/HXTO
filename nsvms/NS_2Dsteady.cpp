#include "NS_2Dsteady.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

typedef unsigned int uint;

double MatrixDet(double dxdt[2][2])
{
	double det = dxdt[0][0] * dxdt[1][1] - dxdt[0][1] * dxdt[1][0];
	
	return det;
}

void Matrix2DInverse(double dxdt[2][2], double dtdx[2][2])
{
	double det = MatrixDet(dxdt);
	dtdx[0][0] = 1.0 / det*(dxdt[1][1]);
	dtdx[0][1] = 1.0 / det*(-dxdt[0][1]);
	dtdx[1][0] = 1.0 / det*(-dxdt[1][0]);
	dtdx[1][1] = 1.0 / det*(dxdt[0][0]);
}

NS_2Dsteady::NS_2Dsteady()
{
	comm = MPI_COMM_WORLD;
	mpiErr = MPI_Comm_rank(comm, &comRank);
	mpiErr = MPI_Comm_size(comm, &comSize);
	nProcess = comSize;	
}

void NS_2Dsteady::BodyForce(double x, double y, double &Fx, double &Fy)
{
	Fx = 0;
	Fy = 0;
}

void NS_2Dsteady::ReadBezierElementProcess(string fn)
{
	string stmp;
	int npts, neles, nfunctions, itmp, itmp1;
	int add(0);

	string fname_cmat = fn + "cmat.txt";
	ifstream fin_cmat;
	fin_cmat.open(fname_cmat);
	if (fin_cmat.is_open())	{
		fin_cmat >> neles;
		//cout << neles << endl;
		// cout << "test 666!!!"<< endl;
		bzmesh_process.resize(ele_process.size());
		for (int i = 0; i<neles; i++){
			if (i == ele_process[add]){
				fin_cmat >> itmp >> nfunctions >> bzmesh_process[add].type;
				bzmesh_process[add].cmat.resize(nfunctions);
				bzmesh_process[add].IEN.resize(nfunctions);
				for (int j = 0; j < nfunctions; j++)
					fin_cmat >> bzmesh_process[add].IEN[j];
				for (int j = 0; j < nfunctions; j++){
					for (int k = 0; k < 16; k++){
						fin_cmat >> bzmesh_process[add].cmat[j][k];
					}
				}
				add++;
			}
			else{
				fin_cmat >> stmp >> nfunctions >> itmp;
				for (int j = 0; j < nfunctions; j++)
					fin_cmat >> stmp;
				for (int j = 0; j < nfunctions; j++)
					for (int k = 0; k < 16; k++)
						fin_cmat >> stmp;
			}
		}
		fin_cmat.close();
		//PetscPrintf(PETSC_COMM_WORLD, "Bezier Matrices Loaded!\n");
	}
	else{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_cmat.c_str());
	}

	string fname_bzpt = fn + "bzpt.txt";
	ifstream fin_bzpt;
	fin_bzpt.open(fname_bzpt);
	add = 0;
	double tmp_xuan;
	if (fin_bzpt.is_open()){
		fin_bzpt >> npts;
		getline(fin_bzpt, stmp);
		for (int e = 0; e < neles; e++)	{
			if (e == ele_process[add]){
				bzmesh_process[add].pts.resize(bzpt_num);
				for (int i = 0; i < bzpt_num; i++){
					fin_bzpt >> bzmesh_process[add].pts[i][0] >> bzmesh_process[add].pts[i][1] >> tmp_xuan; // bzmesh_process[add].pts[i][2];
				}
				add++;
			}
			else{
				for (int i = 0; i < bzpt_num; i++)
					fin_bzpt >> stmp >> stmp >> stmp;
			}
		}
		fin_bzpt.close();
		//PetscPrintf(PETSC_COMM_WORLD, "Bezier Points Loaded!\n");
	}
	else{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname_bzpt.c_str());
	}
}

void NS_2Dsteady::GaussInfo(int ng)
{
	Gpt.clear();
	wght.clear();
	switch (ng)
	{
	case 2:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.2113248654051871;			Gpt[1] = 0.7886751345948129;
		wght[0] = 1.;							wght[1] = 1.;
		break;
	}
	case 3:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.1127016653792583;			Gpt[1] = 0.5;							Gpt[2] = 0.8872983346207417;
		wght[0] = 0.5555555555555556;			wght[1] = 0.8888888888888889;			wght[2] = 0.5555555555555556;
		break;
	}
	case 4:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.06943184420297371;			Gpt[1] = 0.33000947820757187;			Gpt[2] = 0.6699905217924281;			Gpt[3] = 0.9305681557970262;
		wght[0] = 0.3478548451374539;			wght[1] = 0.6521451548625461;			wght[2] = 0.6521451548625461;			wght[3] = 0.3478548451374539;
		break;
	}
	case 5:
	{
		Gpt.resize(ng);
		wght.resize(ng);
		Gpt[0] = 0.046910077030668;			Gpt[1] = 0.2307653449471585;		Gpt[2] = 0.5;						Gpt[3] = 0.7692346550528415;		Gpt[4] = 0.953089922969332;
		wght[0] = 0.2369268850561891;		wght[1] = 0.4786286704993665;		wght[2] = 0.5688888888888889;		wght[3] = 0.4786286704993665;		wght[4] = 0.2369268850561891;
		break;
	}
	default:
	{
		Gpt.resize(2);
		wght.resize(2);
		Gpt[0] = 0.2113248654051871;			Gpt[1] = 0.7886751345948129;
		wght[0] = 1.;							wght[1] = 1.;
		break;
	}
	}
}

void NS_2Dsteady::InitializeProblem(const int ndof, const int n_bz, const vector<double>& Vel0, const vector<double>& Pre0, const vector<double>& var)
{
	
	MPI_Barrier(comm);

	PetscPrintf(PETSC_COMM_WORLD, "Initializing...\n");
	
	/*Initialize parameters*/
	GaussInfo(4);
	n_bzmesh = n_bz;
	Vel.resize(ndof * 2);
	Pre.resize(ndof);
	Vel = Vel0;
	Pre = Pre0;
	if (var.size() != 0)
	{
		nu = 0.25;
		rou = 0.5;
		alphaM = 0.5*(3 - rou) / (1 + rou);
		alphaF = 1 / (1 + rou);
		Gama = 0.5 + alphaM - alphaF;
		velocity_max = var[1];
	}
	else
	{
		cerr << "0 variables!\n"; getchar();
	}
	fx = 0;
	fy = 0;
	// fz = 0;

	/*Initialize petsc vector, matrix*/
	PetscInt mat_dim = ndof * 3;
	ierr = MatCreate(PETSC_COMM_WORLD, &GK); 
	ierr = MatSetSizes(GK, PETSC_DECIDE, PETSC_DECIDE, mat_dim, mat_dim);
	ierr = MatSetType(GK, MATMPIAIJ);
	ierr = MatMPIAIJSetPreallocation(GK, 500, NULL, 500, NULL);
	MatGetOwnershipRange(GK, &rstart, &rend);		
	ierr = MatSetOption(GK, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(GK); 
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &GR);
	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, mat_dim, &temp_solution);
}

void NS_2Dsteady::BasisFunction(double u, double v, int nen, const vector<array<double, 2>>& pt, const vector<array<double, 16>> &cmat, vector<double> &Nx, vector<array<double, 2>> &dNdx, vector<array<array<double, 2>, 2>> &dN2dx2, double dudx[2][2], double& detJ)
{
	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
	// double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };

	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
	// double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
	
	double dN2du2[4] = { 6.*(1. - u),-12. + 18.*u,6. - 18.*u,6.*u };
	double dN2dv2[4] = { 6.*(1. - v),-12. + 18.*v,6. - 18.*v,6.*v };
	// double dN2dw2[4] = { 6.*(1. - w),-12. + 18.*w,6. - 18.*w,6.*w };

	double dNdt[bzpt_num][2];
	double dN2dt2[bzpt_num][2][2];
	double Nx_bz[bzpt_num];
	double dNdx_bz[bzpt_num][2];
	double dN2dx2_bz[bzpt_num][2][2];

	Nx.clear();
	dNdx.clear();
	dN2dx2.clear();
	Nx.resize(nen, 0);
	dNdx.resize(nen, { 0 });
	dN2dx2.resize(nen, { {0} });

	int i, j, k, a, b, c, loc;
	loc = 0;
	for (j = 0; j < 4; j++)
	{
		for (k = 0; k < 4; k++)
		{
			Nx_bz[loc] = Nu[k] * Nv[j];
			dNdt[loc][0] = dNdu[k] * Nv[j];
			dNdt[loc][1] = Nu[k] * dNdv[j];
			dN2dt2[loc][0][0] = dN2du2[k] * Nv[j];
			dN2dt2[loc][0][1] = dNdu[k] * dNdv[j];
			dN2dt2[loc][1][0] = dNdu[k] * dNdv[j];
			dN2dt2[loc][1][1] = Nu[k] * dN2dv2[j];
			loc++;
		}
	}

	double dxdt[2][2] = { {0} };
	for (loc = 0; loc < bzpt_num; loc++)
	{
		for (a = 0; a<2; a++)
		{
			for (b = 0; b<2; b++)
			{
 				dxdt[a][b] += pt[loc][a] * dNdt[loc][b];
			}
		}	
			
	}

	double dtdx[2][2] = { { 0 } };
	Matrix2DInverse(dxdt, dtdx);

	//1st derivatives
	for (i = 0; i<bzpt_num; i++)
	{
		dNdx_bz[i][0] = dNdt[i][0] * dtdx[0][0] + dNdt[i][1] * dtdx[1][0];
		dNdx_bz[i][1] = dNdt[i][0] * dtdx[0][1] + dNdt[i][1] * dtdx[1][1];
	}
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			dudx[i][j] = dtdx[i][j];

	detJ = MatrixDet(dxdt);
	detJ = 0.25*detJ;

	//2nd derivatives
	double dx2dt2[2][4] = { {0} };
	double dt2dx2[2][4] = { {0} };
	for (int l = 0; l < 2; l++)	
		for (loc = 0; loc<bzpt_num; loc++)
			for (a = 0; a<2; a++)
				for (b = 0; b<2; b++)
					dx2dt2[l][2*a+b] += pt[loc][l] * dN2dt2[loc][a][b];

	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			for (k = 0; k < 2; k++)
				for (a = 0; a < 2; a++)
					for (b = 0; b < 2; b++)
						for (c = 0; c < 2; c++)
							dt2dx2[c][2*i+j] -= dx2dt2[k][2*a+b]*dtdx[a][i]*dtdx[b][j]*dtdx[c][k];


	
	for (loc = 0; loc < bzpt_num; loc++)
		for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++)
				dN2dx2_bz[loc][i][j] = 0.;

	for (loc = 0; loc < bzpt_num; loc++){
		for (i = 0; i < 2; i++)	{
			for (j = 0; j < 2; j++)	{
				for (a = 0; a<2; a++){
					for (b = 0; b<2; b++){
						dN2dx2_bz[loc][i][j] += dN2dt2[loc][a][b] * dtdx[a][i]*dtdx[b][j];
					}
					dN2dx2_bz[loc][i][j] += dNdt[loc][a] * dt2dx2[a][2*i+j];
				}
			}
		}
	}
	
	for (i = 0; i < nen; i++){
		for (j = 0; j < bzpt_num; j++)	{
			Nx[i] += cmat[i][j] * Nx_bz[j];
			for (int m = 0; m < 2; m++)	{
				dNdx[i][m] += cmat[i][j] * dNdx_bz[j][m];
				for (int n = 0; n < 2; n++)	{
					dN2dx2[i][m][n] += cmat[i][j] * dN2dx2_bz[j][m][n];
				}
			}
		}
	}

}

void NS_2Dsteady::PointFormValue(vector<double> &Nx, const vector<array<double, 3>> &U, double Value[3])
{
	for (int i = 0; i < 3; i++)
		Value[i] = 0;

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < Nx.size(); j++)
			Value[i] += U[j][i] * Nx[j];
}

void NS_2Dsteady::PointFormGrad(vector<array<double, 2>>& dNdx, const vector<array<double, 3>> &U, double Value[3][2])
{
	for (int i = 0; i < 3; i++)	
		for (int j = 0; j < 2; j++)	
			Value[i][j] = 0.;

	for (int i = 0; i < 3; i++) 
		for (int j = 0; j < 2; j++)	
			for (int k = 0; k < dNdx.size(); k++) 
				Value[i][j] += U[k][i] * dNdx[k][j];
}

void NS_2Dsteady::PointFormHess(vector<array<array<double, 2>, 2>>& d2Ndx2, const vector<array<double, 3>> &U, double Value[3][2][2])
{
	for (int i = 0; i < 3; i++)	{
		for (int j = 0; j < 2; j++)	{
			for (int k = 0; k < 2; k++)	{
				Value[i][j][k] = 0.;
			}
		}
	}
	for (int i = 0; i < 3; i++)	{
		for (int j = 0; j < 2; j++)	{
			for (int k = 0; k < 2; k++)	{
				for (int l = 0; l < d2Ndx2.size(); l++)	{
					Value[i][j][k] += U[l][i] * d2Ndx2[l][j][k];
				}
			}
		}
	}
}

void NS_2Dsteady::Tau(double J[2][2], double u[3], double &tauM, double &tauC) 
{
	/*calculate stabilization parameter*/
	double C_I = 1.0 / 12.0;

	int i, j, k;

	double G[2][2] = { { 0 } };
	for (i = 0; i<2; i++)
		for (j = 0; j<2; j++)
			for (k = 0; k<2; k++)
				G[i][j] += J[k][i] * J[k][j];

	double g[2] = { 0 };
	for (i = 0; i<2; i++)
		for (j = 0; j<2; j++)
			g[i] += J[j][i];

	double G_G = 0;
	for (i = 0; i<2; i++)
		for (j = 0; j<2; j++)
			G_G += G[i][j] * G[i][j];

	double g_g = 0;
	for (i = 0; i<2; i++)
		g_g += g[i] * g[i];

	double u_G_u = 0;
	for (i = 0; i<2; i++)
		for (j = 0; j<2; j++)
			u_G_u += u[i] * G[i][j] * u[j];

	
	tauM = u_G_u + C_I * nu*nu * G_G;
	tauM = 1 / sqrt(tauM);
	
	tauC = (tauM)* g_g;
	tauC = 1 / (tauC);
}

void NS_2Dsteady::FineScale(double tauM, double tauC, double u[2], double u_x[2], double u_y[2], double u_xx[2], double u_yy[2], double p, double p_x, double p_y, double u_s[2], double &p_s)
{
	/*calculate fine scale in VMS*/
	u_s[0] = (u[0] * u_x[0] + u[1] * u_y[0]) + p_x - nu*(u_xx[0] + u_yy[0]) - fx;
	u_s[1] = (u[0] * u_x[1] + u[1] * u_y[1]) + p_y - nu*(u_xx[1] + u_yy[1]) - fy;
	// u_s[0] = (u[0] * u_x[0] + u[1] * u_y[0] + u[2] * u_z[0]) + p_x - nu*(u_xx[0] + u_yy[0] + u_zz[0]) - fx;
	// u_s[1] = (u[0] * u_x[1] + u[1] * u_y[1] + u[2] * u_z[1]) + p_y - nu*(u_xx[1] + u_yy[1] + u_zz[1]) - fy;
	// u_s[2] = (u[0] * u_x[2] + u[1] * u_y[2] + u[2] * u_z[2]) + p_z - nu*(u_xx[2] + u_yy[2] + u_zz[2]) - fz;
	
	p_s = u_x[0] + u_y[1];
	// p_s = u_x[0] + u_y[1] + u_z[2];
	
	u_s[0] *= -tauM;
	u_s[1] *= -tauM;
	// u_s[2] *= -tauM;
	p_s *= -tauC;
}

void NS_2Dsteady::Residual(vector<double>& Nx, vector<array<double, 2>>& dNdx, vector<array<array<double, 2>, 2>>& dN2dx2, double dudx[2][2], const double detJ, const vector<array<double, 3>> &U, vector<array<double, 3>> Re)
{
	/*calculate residual of the equation*/
	double U_t[3], UU[3];
	double grad_U[3][2];
	double der2_U[3][2][2];

	PointFormValue(Nx, U, UU);
	PointFormGrad(dNdx, U, grad_U);
	PointFormHess(dN2dx2, U, der2_U);

	double u[2] = { UU[0],UU[1] };
	double u_x[2] = { grad_U[0][0],grad_U[1][0] };
	double u_y[2] = { grad_U[0][1],grad_U[1][1] };
	// double u_z[3] = { grad_U[0][2],grad_U[1][2],grad_U[2][2] };
	double u_xx[2] = { der2_U[0][0][0],der2_U[1][0][0] };
	double u_yy[2] = { der2_U[0][1][1],der2_U[1][1][1] };
	// double u_zz[3] = { der2_U[0][2][2],der2_U[1][2][2],der2_U[2][2][2] };
	double p = UU[2];
	double p_x = grad_U[2][0], p_y = grad_U[2][1];


	double InvGradMap[2][2];
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			InvGradMap[i][j] = dudx[i][j];
		}
	}

	double tauM, tauC;
	Tau(InvGradMap, UU, tauM, tauC);
	double u_s[2], p_s;
	FineScale(tauM, tauC, u, u_x, u_y, u_xx, u_yy, p, p_x, p_y, u_s, p_s);

	int a, nen = Nx.size();
	for (a = 0; a<nen; a++) {
		double Na = Nx[a];
		double Na_x = dNdx[a][0];
		double Na_y = dNdx[a][1];
		/* ----- */
		double Rux, Ruy, Rp;
		// -L(W)
		Rux = -Na*fx;
		Ruy = -Na*fy;
		Rp = 0.0;
		// += B_1(W,U)
		Rux += -Na_x*p + nu*(Na_x*(u_x[0] + u_x[0]) + Na_y*(u_y[0] + u_x[1]));
		Ruy += -Na_y*p + nu*(Na_x*(u_x[1] + u_y[0]) + Na_y*(u_y[1] + u_y[1]));
		Rp  += +Na*(u_x[0] + u_y[1]);
		// += Btilde_1(W,U')
		Rux += -(Na_x*p_s);
		Ruy += -(Na_y*p_s);
		Rp  += -(Na_x*u_s[0] + Na_y*u_s[1]);
		// += B_2(W,U,U+U') [after integration by parts, div(u+u')=0]
		Rux += +Na * ((u[0] + u_s[0])*u_x[0] + (u[1] + u_s[1])*u_y[0]);
		Ruy += +Na * ((u[0] + u_s[0])*u_x[1] + (u[1] + u_s[1])*u_y[1]);
		// += B_2(W,U',U+U')
		Rux += -(Na_x*u_s[0] * (u[0] + u_s[0]) + Na_y*u_s[0] * (u[1] + u_s[1]));
		Ruy += -(Na_x*u_s[1] * (u[0] + u_s[0]) + Na_y*u_s[1] * (u[1] + u_s[1]));
		/* ----- */
		Re[a][0] += Rux * detJ;
		Re[a][1] += Ruy * detJ;
		Re[a][2] += Rp  * detJ;
	}	
	/* Following is for 3D case
	for (a = 0; a<nen; a++) {
		double Na = Nx[a];
		double Na_x = dNdx[a][0];
		double Na_y = dNdx[a][1];
		double Na_z = dNdx[a][2];
		
		double Rux, Ruy, Ruz, Rp;
		
		Rux = -Na*fx;
		Ruy = -Na*fy;
		Ruz = -Na*fz;
		Rp = 0.0;
		
		Rux += - Na_x*p + nu*(Na_x*(u_x[0] + u_x[0]) + Na_y*(u_y[0] + u_x[1]) + Na_z*(u_z[0] + u_x[2]));
		Ruy += - Na_y*p + nu*(Na_x*(u_x[1] + u_y[0]) + Na_y*(u_y[1] + u_y[1]) + Na_z*(u_z[1] + u_y[2]));
		Ruz += - Na_z*p + nu*(Na_x*(u_x[2] + u_z[0]) + Na_y*(u_y[2] + u_z[1]) + Na_z*(u_z[2] + u_z[2]));
		Rp += Na*(u_x[0] + u_y[1] + u_z[2]);
		
		Rux += -(Na_x*p_s);
		Ruy += -(Na_y*p_s);
		Ruz += -(Na_z*p_s);
		Rp += -(Na_x*u_s[0] + Na_y*u_s[1] + Na_z*u_s[2]);
		
		Rux += +Na * ((u[0] + u_s[0])*u_x[0] + (u[1] + u_s[1])*u_y[0] + (u[2] + u_s[2])*u_z[0]);
		Ruy += +Na * ((u[0] + u_s[0])*u_x[1] + (u[1] + u_s[1])*u_y[1] + (u[2] + u_s[2])*u_z[1]);
		Ruz += +Na * ((u[0] + u_s[0])*u_x[2] + (u[1] + u_s[1])*u_y[2] + (u[2] + u_s[2])*u_z[2]);
		
		Rux += -(Na_x*u_s[0] * (u[0] + u_s[0]) + Na_y*u_s[0] * (u[1] + u_s[1]) + Na_z*u_s[0] * (u[2] + u_s[2]));
		Ruy += -(Na_x*u_s[1] * (u[0] + u_s[0]) + Na_y*u_s[1] * (u[1] + u_s[1]) + Na_z*u_s[1] * (u[2] + u_s[2]));
		Ruz += -(Na_x*u_s[2] * (u[0] + u_s[0]) + Na_y*u_s[2] * (u[1] + u_s[1]) + Na_z*u_s[2] * (u[2] + u_s[2]));
		
		Re[a][0] += Rux * detJ;
		Re[a][1] += Ruy * detJ;
		Re[a][2] += Ruz * detJ;
		Re[a][3] += Rp  * detJ;
	}
	*/
}

void NS_2Dsteady::Tangent(vector<double> &Nx, vector<array<double, 2>>& dNdx, double dudx[2][2], const double detJ, const vector<array<double, 3>>& U, vector<array<vector<array<double, 3>>, 3>>& Ke)
{
	/*calculate tangent matrix*/
	double u[3];
	PointFormValue(Nx, U, u);
	double ux = u[0];
	double uy = u[1];
	// double uz = u[2];

	double InvGradMap[2][2];
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			InvGradMap[i][j] = dudx[i][j];
		}
	}

	double tauM, tauC;
	Tau(InvGradMap, u, tauM, tauC);
	int a, b, nen = Nx.size();
	for (a = 0; a<nen; a++) {
		double Na = Nx[a];
		double Na_x = dNdx[a][0];
		double Na_y = dNdx[a][1];
		// double Na_z = dNdx[a][2];
		for (b = 0; b<nen; b++) {
			double Nb = Nx[b];
			double Nb_x = dNdx[b][0];
			double Nb_y = dNdx[b][1];
			// double Nb_z = dNdx[b][2];
			/* ----- */
			int i, j;
			double T[3][3];
			double Tii =
				(
					+ Na * (ux * Nb_x + uy * Nb_y)
					+ nu * (Na_x * Nb_x + Na_y * Nb_y)
					+ tauM * (ux * Na_x + uy * Na_y) *
					/**/     ((ux * Nb_x + uy * Nb_y))
					);
			
			T[0][0] = (+nu * Na_x * Nb_x + tauC * Na_x * Nb_x);
			T[0][1] = (+nu * Na_y * Nb_x + tauC * Na_x * Nb_y);
			// T[0][2] = (+nu * Na_z * Nb_x + tauC * Na_x * Nb_z);
			//			
			T[1][0] = (+nu * Na_x * Nb_y + tauC * Na_y * Nb_x);
			T[1][1] = (+nu * Na_y * Nb_y + tauC * Na_y * Nb_y);
			// T[1][2] = (+nu * Na_z * Nb_y + tauC * Na_y * Nb_z);
			//			
			// T[2][0] = (+nu * Na_x * Nb_z + tauC * Na_z * Nb_x);
			// T[2][1] = (+nu * Na_y * Nb_z + tauC * Na_z * Nb_y);
			// T[2][2] = (+nu * Na_z * Nb_z + tauC * Na_z * Nb_z);
			T[0][0] += Tii;
			T[1][1] += Tii;
			// T[2][2] += Tii;
			
			T[0][2] = -Na_x * Nb + tauM * (ux * Na_x + uy * Na_y) * Nb_x;
			T[1][2] = -Na_y * Nb + tauM * (ux * Na_x + uy * Na_y) * Nb_y;
			// T[2][3] = -Na_z * Nb + tauM * (ux * Na_x + uy * Na_y + uz * Na_z) * Nb_z;
			
			T[2][0] = +Na * Nb_x + tauM * Na_x * ((ux * Nb_x + uy * Nb_y));
			T[2][1] = +Na * Nb_y + tauM * Na_y * ((ux * Nb_x + uy * Nb_y));
			// T[3][2] = +Na * Nb_z + tauM * Na_z * ((ux * Nb_x + uy * Nb_y + uz * Nb_z));
			
			T[2][2] = +tauM * (Na_x * Nb_x + Na_y * Nb_y);
			
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					Ke[a][i][b][j] += T[i][j] * detJ;
		}
	}
}

void NS_2Dsteady::MatrixAssembly(vector<array<vector<array<double, 3>>, 3>>& Ke, const vector<int>& IEN, Mat& GK)
{
	int i, j, A, B, m, n;	
	int row_start, row_end, row_now;
	int add=0;

	PetscInt *nodeList = new PetscInt[IEN.size()*3];
	PetscReal *tmpGK = new PetscReal[IEN.size() * 3 * IEN.size() * 3];
	
	for (m = 0; m<IEN.size(); m++){
		A = IEN[m];
		for (i = 0; i < 3; i++){
			nodeList[3 * m + i] = 3 * A + i;
			for (n = 0; n<IEN.size(); n++){
				B = IEN[n];
				for (j = 0; j < 3; j++){
					tmpGK[add] = Ke[m][i][n][j];
					add++;								 
				}
			}
		}
	}		
	MatSetValues(GK, IEN.size() * 3, nodeList, IEN.size() * 3, nodeList, tmpGK, ADD_VALUES);
	delete nodeList;
	delete tmpGK;	
}

void NS_2Dsteady::ResidualAssembly(vector<array<double, 3>> &Re, const vector<int>& IEN, Vec& GR)
{
	int i, j, A, B, m, n;
	int add = 0;


	PetscInt *nodeList = new PetscInt[IEN.size() * 3];
	PetscReal *tmpGR = new PetscReal[IEN.size() * 3];
	for (m = 0; m<IEN.size(); m++)
	{
		A = IEN[m];
		for (i = 0; i < 3; i++)	{
			nodeList[3 * m + i] = 3 * A + i;
			tmpGR[add] = -Re[m][i];
			add++;
		}
	}
	
	VecSetValues(GR, IEN.size() * 3, nodeList, tmpGR, ADD_VALUES);
	delete  nodeList;
	delete  tmpGR;
}

void NS_2Dsteady::BuildLinearSystemProcess(const vector<Vertex2D>& cpts, const vector<array<double, 2>>& velocity_bc, const vector<double> velocity_node, const vector<double> pressure_node)
{
	/*Build linear system in each process*/
	int e;
	//cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for "<< bzmesh_process.size() <<" elements.\n";
	for (e=0;e<bzmesh_process.size();e++){
		int nen, A;
		double ux_bc, uy_bc, p_bc;
	
		double dudx[2][2];
		double detJ;
		vector<double> Nx;
		vector<array<double, 2>> dNdx;
		vector<array<array<double, 2>, 2>> dN2dx2;
		vector<array<double, 3>> Re;
		vector<array<vector<array<double, 3>>, 3>> Ke;
	
		vector<array<double, 3>> v_node;		
		
		nen = bzmesh_process[e].IEN.size();
		
		Nx.clear(); Nx.resize(nen, 0);
		dNdx.clear(); dNdx.resize(nen, { 0 });
		dN2dx2.clear(); dN2dx2.resize(nen, { {0} });
		Re.clear(); Re.resize(nen);
		Ke.clear(); Ke.resize(nen);
		v_node.clear(); v_node.resize(nen);
	
		for (int i = 0; i < nen; i++)
			for (int j = 0; j < 3; j++)
				Ke[i][j].resize(nen);
	
		for (int i = 0; i < nen; i++){
			for (int m = 0; m < 3; m++)	{
				for (int j = 0; j < nen; j++){
					for (int n = 0; n < 3; n++)	{
						Ke[i][m][j][n] = 0.;
					}	
				}
				Re[i][m] = 0.;
			}
		}
	
		for (int i = 0; i < nen; i++){
			for (int j = 0; j < 2; j++){
				v_node[i][j] = velocity_max * velocity_node[bzmesh_process[e].IEN[i] * 2 + j];
			}
			v_node[i][2] = pressure_node[bzmesh_process[e].IEN[i]];
		}
	
		for (int i = 0; i < Gpt.size(); i++)
		{
			for (int j = 0; j < Gpt.size(); j++)
			{
				BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
				detJ = wght[i] * wght[j] * detJ;
				Tangent(Nx, dNdx, dudx, detJ, v_node, Ke);
				Residual(Nx, dNdx, dN2dx2, dudx, detJ, v_node, Re);
			}
		}
	
		/*Apply Boundary Condition*/
		for (int i = 0; i < nen; i++){
			A = bzmesh_process[e].IEN[i];
			if (cpts[A].label == 1) {
				//inlet
				ux_bc = velocity_max*velocity_bc[A][0] - velocity_node[A * 2];
				// cout << "Xuan_Test_666!!!" << ux_bc << endl;				
				ApplyBoundaryCondition(ux_bc, i, 0, Ke, Re);
				uy_bc = velocity_max*velocity_bc[A][1] - velocity_node[A * 2 + 1];
				// cout << "Xuan_Test_999!!!"<< uy_bc << endl;
				ApplyBoundaryCondition(uy_bc, i, 1, Ke, Re);
				// uz_bc = velocity_max*velocity_bc[A][2] - velocity_node[A * 3 + 2];
				// ApplyBoundaryCondition(uz_bc, i, 2, Ke, Re);
			}
			if (cpts[A].label == 0) {
				//wall
				ux_bc = 0.0 - velocity_node[A * 2];
				ApplyBoundaryCondition(ux_bc, i, 0, Ke, Re);
        // cout << "Xuan_Test_999!!!"<< ux_bc << endl;
				uy_bc = 0.0 - velocity_node[A * 2 + 1];
				ApplyBoundaryCondition(uy_bc, i, 1, Ke, Re);
				// uz_bc = 0.0 - velocity_node[A * 3 + 2];
				// ApplyBoundaryCondition(uz_bc, i, 2, Ke, Re);
			}
			if (cpts[A].label == 3) {
				//wall
				ux_bc = 0.0 - velocity_node[A * 2];
				ApplyBoundaryCondition(ux_bc, i, 0, Ke, Re);
				uy_bc = 0.0 - velocity_node[A * 2 + 1];
				ApplyBoundaryCondition(uy_bc, i, 1, Ke, Re);
				// uz_bc = 0.0 - velocity_node[A * 3 + 2];
				// ApplyBoundaryCondition(uz_bc, i, 2, Ke, Re);
			}		
		}

		/*Start element vector assembly*/
		ResidualAssembly(Re, bzmesh_process[e].IEN, GR);
	
		/*Start element matrix assembly*/
		MatrixAssembly(Ke, bzmesh_process[e].IEN, GK);
	}
	//cout << "Process " << comRank << " :complete build matrix and vector!\n";
	VecAssemblyBegin(GR);
	MatAssemblyBegin(GK, MAT_FINAL_ASSEMBLY);	
}

void NS_2Dsteady::ApplyBoundaryCondition(const double bc_value, int pt_num, int variable_num, vector<array<vector<array<double, 3>>, 3>>& Ke, vector<array<double, 3>> &Re)
{
	int j, k;
	for (j = 0; j < Re.size(); j++)	{
		for (k = 0; k < 3; k++)		{
			Re[j][k] += bc_value*Ke[j][k][pt_num][variable_num];
		}
	}
	for (j = 0; j < Re.size(); j++)	{
		for (k = 0; k < 3; k++)		{
			Ke[pt_num][variable_num][j][k] = 0.0;		Ke[j][k][pt_num][variable_num] = 0.0;
		}
	}
	Re[pt_num][variable_num] = -bc_value;
	Ke[pt_num][variable_num][pt_num][variable_num] = 1.0;
}

void NS_2Dsteady::AssignProcessor(vector<vector<int>> &ele_proc)
{
	/*Assign the partitioned bezier elements to this processor */
	for (int i = 0; i < ele_proc[comRank].size(); i++)
		ele_process.push_back(ele_proc[comRank][i]);
}

void DebugVisualizeMat(Mat mat_debug, string fname){
	PetscErrorCode ierr;
	PetscViewer viewer;
	PetscViewerFormat format = PETSC_VIEWER_ASCII_MATLAB;

	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
	ierr = PetscViewerPushFormat(viewer, format);
	ierr = MatView(mat_debug, viewer);
	ierr = PetscViewerPopFormat(viewer);
	ierr = PetscViewerDestroy(&viewer);
}

void DebugVisualizeVec(Vec vec_debug, string fname){
	PetscErrorCode ierr;
	PetscViewer viewer;
	PetscViewerFormat format = PETSC_VIEWER_ASCII_MATLAB;

	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname.c_str(), &viewer);
	ierr = PetscViewerPushFormat(viewer, format);
	ierr = VecView(vec_debug, viewer);
	ierr = PetscViewerPopFormat(viewer);
	ierr = PetscViewerDestroy(&viewer);
}

void NS_2Dsteady::Pressure_Drop(const vector<Vertex2D>& cpts, const vector<double> pressure_node, double& pre_drop)
{
	/*Compute pressure drop between outlet and inlet*/
	int e;
	//cout << "Process:" << comRank << " out of " << nProcess << " Start Loop for "<< bzmesh_process.size() <<" elements.\n";
	double pd_inlet, pd_outlet, pd_tmp = 0.0;
  	double area_inlet, area_outlet;
  	double area = 0.0155156;
	for (e=0;e<bzmesh_process.size();e++){
		int nen, A;
		// double ux_bc, uy_bc, p_bc;
	
		double dudx[2][2];
    	double dxdt[2][2];
		double detJ;
		double det_bc_J1, det_bc_J2;
		vector<double> Nx;
		vector<array<double, 2>> dNdx;
		vector<array<array<double, 2>, 2>> dN2dx2;
	
		vector<array<double, 3>> v_node;		
		
		nen = bzmesh_process[e].IEN.size();
		
		Nx.clear(); Nx.resize(nen, 0);
		dNdx.clear(); dNdx.resize(nen, { 0 });
		dN2dx2.clear(); dN2dx2.resize(nen, { {0} });

		/*Check boundary flags*/
		int flag_bc = 0;
		for (int i = 0; i < nen; i++){
			int A = bzmesh_process[e].IEN[i];
			if (cpts[A].label == 1)	{
				// inlet
				// N0_bc = par[9];
				flag_bc = 1;
				break;
			}
			if (cpts[A].label == 2)	{
				//outlet
				// N0_bc = par[10];
				flag_bc = 2;
				break;
				// Nplus_bc = par[10];
			}
		}
		/*Boundary integral for pressure integration*/
		if(flag_bc == 1){
			for (int i = 0; i < Gpt.size(); i++)
			{
				for (int j = 0; j < Gpt.size(); j++)
				{
					BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
					//Matrix2DInverse(dudx, dxdt);
					/*
					for(int k = 0; k < nen; k++){
						int A = bzmesh_process[e].IEN[k];
						if(cpts[A].label == 1){
							det_bc_J1 += dNdx[k][1]*dNdx[k][1];
						}
					} */
					//det_bc_J1 = sqrt(dxdt[0][1]*dxdt[0][1] + dxdt[1][1]*dxdt[1][1]);
					// det_bc_J1 = wght[i] * wght[j] * sqrt(det_bc_J1);
					//Tangent(Nx, dNdx, dudx, detJ, v_node, Ke);
					//Residual(Nx, dNdx, dN2dx2, dudx, detJ, v_node, Re);
          			detJ = wght[i] * wght[j] * detJ;
					int n;
					for(n = 0; n < nen; n++){
						int A = bzmesh_process[e].IEN[n];
						if(flag_bc == 1){
							pd_inlet += pressure_node[bzmesh_process[e].IEN[n]] * Nx[n] *detJ;						
						}
					}
              		area_inlet += detJ;	
  			      	det_bc_J1 = 0.0;  
				}
			}
      		flag_bc = 0;
		}
		if(flag_bc == 2){
			for (int i = 0; i < Gpt.size(); i++)
			{
				for (int j = 0; j < Gpt.size(); j++)
				{
					BasisFunction(Gpt[i], Gpt[j], nen, bzmesh_process[e].pts, bzmesh_process[e].cmat, Nx, dNdx, dN2dx2, dudx, detJ);
					/*
					for(int k = 0; k < nen; k++){
						int A = bzmesh_process[e].IEN[k];
						if(cpts[A].label == 2){
                			det_bc_J2 += dNdx[k][1]*dNdx[k][1];
					 	}
					} */
					// Matrix2DInverse(dudx, dxdt);
					// det_bc_J2 = sqrt(dxdt[0][1]*dxdt[0][1] + dxdt[1][1]*dxdt[1][1]) * 0.5;
					// det_bc_J2 = wght[i] * wght[j] * sqrt(det_bc_J2);
					//Tangent(Nx, dNdx, dudx, detJ, v_node, Ke);
					//Residual(Nx, dNdx, dN2dx2, dudx, detJ, v_node, Re);
          			detJ = wght[i] * wght[j] * detJ;
					int n;
					for(n = 0; n < nen; n++){
						int A = bzmesh_process[e].IEN[n];
						if(flag_bc == 2){
							pd_outlet += pressure_node[bzmesh_process[e].IEN[n]] * Nx[n] * detJ;						
						}
					}
          			area_outlet += detJ;	
					det_bc_J2 = 0.0;	
				}
			}
      		flag_bc = 0;
		}
	}
	/*Compute the pressure drop*/
	pd_tmp = (pd_inlet - pd_outlet)/area;
	pre_drop = pd_tmp;
	//cout << "test_xuan: Pressure Drop Calculation:   " << pre_drop << endl;
}

void NS_2Dsteady::Run(const vector<Vertex2D>& cpts, const vector<Element2D>& tmesh, const vector<array<double, 2>>& velocity_bc, string fn)
{
	
	int n_iterator(1);
	int l(0);
	time_t t0, t1;	
	vector<double> V_delta(2 * cpts.size()), P_delta(cpts.size());
	
	ReadBezierElementProcess(fn);

	for (l = 0; l<n_iterator; l++){
	
		//PetscPrintf(PETSC_COMM_WORLD, "Iteration Step: %d\n", l);
		/*Build Linear System*/
		//PetscPrintf(PETSC_COMM_WORLD, "Building Linear System...\n");	
		MatZeroEntries(GK);
		VecSet(GR, 0);		
		t0 = time(NULL);
		BuildLinearSystemProcess(cpts, velocity_bc, Vel, Pre);
		VecAssemblyEnd(GR);
		//PetscPrintf(PETSC_COMM_WORLD, "Done Vector Assembly...\n");
		MatAssemblyEnd(GK, MAT_FINAL_ASSEMBLY);	
		t1 = time(NULL);
		//PetscPrintf(PETSC_COMM_WORLD, "Done Matrix Assembly with time: %d \n", t1 - t0);		
		MPI_Barrier(comm);
		PetscPrintf(PETSC_COMM_WORLD, "Solving...\n");

		t0 = time(NULL);
		/*Petsc KSP solver setting*/
		KSPCreate(PETSC_COMM_WORLD, &ksp);
		KSPReset(ksp);
		KSPSetOperators(ksp, GK, GK);		
		KSPGetPC(ksp, &pc);
		PCSetType(pc, PCBJACOBI);
		KSPSetType(ksp, KSPGMRES);
		KSPGMRESSetRestart(ksp, 500);
		KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
		KSP *subksp;
		PC subpc;
		PetscInt first,nlocal;
		KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, 100000);
		KSPSetPCSide(ksp, PC_RIGHT);
		KSPSetFromOptions(ksp);
		KSPSetUp(ksp);
		PCBJacobiGetSubKSP(pc, &nlocal, &first, &subksp);
		for (int i = 0; i<nlocal; i++) {
			KSPGetPC(subksp[i], &subpc);		
			PCSetType(subpc, PCILU);
			KSPSetType(subksp[i], KSPGMRES);
			KSPSetInitialGuessNonzero(subksp[i], PETSC_TRUE);
			KSPSetPCSide(subksp[i], PC_RIGHT);
		}
	
		/*Solving the equation*/
		KSPSolve(ksp, GR, temp_solution);
		//KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
		PetscPrintf(PETSC_COMM_WORLD, "------------------------------\n");
		PetscInt its;
		KSPGetIterationNumber(ksp, &its);
		PetscPrintf(PETSC_COMM_WORLD, "iterations %d\n", its);
		KSPConvergedReason reason;
		KSPGetConvergedReason(ksp, &reason);
		PetscPrintf(PETSC_COMM_WORLD, "KSPConvergedReason: %D\n", reason);
		t1 = time(NULL);
		PetscPrintf(PETSC_COMM_WORLD, "Done Solving with time: %d \n", t1 - t0);

		// DebugVisualizeMat(GK, "./example/GKMat.m");
		// DebugVisualizeVec(GR, "./example/GRMat.m");

		/*Collect the solution from all processors*/
		Vec temp_solution_seq;
		VecScatter ctx;
		PetscScalar    *_a;
		VecScatterCreateToAll(temp_solution, &ctx, &temp_solution_seq);
		VecScatterBegin(ctx, temp_solution, temp_solution_seq, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(ctx, temp_solution, temp_solution_seq, INSERT_VALUES, SCATTER_FORWARD);
		VecGetArray(temp_solution_seq, &_a);		
		MPI_Barrier(comm);
		
		for (uint i = 0; i < cpts.size(); i++){
			V_delta[2 * i] = PetscRealPart(_a[3 * i]);
			V_delta[2 * i + 1] = PetscRealPart(_a[3 * i + 1]);
			// V_delta[3 * i + 2] = PetscRealPart(_a[4 * i + 2]);
			P_delta[i] = PetscRealPart(_a[3 * i + 2]);
		}
		VecRestoreArray(temp_solution_seq, &_a);
		VecScatterDestroy(&ctx);
		VecDestroy(&temp_solution_seq);
		
		for (uint i = 0; i < cpts.size(); i++){
			Vel[2 * i] += V_delta[2 * i];
			Vel[2 * i + 1] += V_delta[2 * i + 1];
			// Vel[3 * i + 2] += V_delta[3 * i + 2];
			Pre[i] += P_delta[i];
		}

		double pd_xuan(0.0);
		Pressure_Drop(cpts, Pre, pd_xuan);
		//cout << "Test_777:   " << pd_xuan << endl;
   
		double * tmp_xuan = NULL;
		if(comRank == 0){
			tmp_xuan = (double*)malloc(sizeof(double)*nProcess);
		}
		MPI_Gather(&pd_xuan, 1, MPI_DOUBLE, tmp_xuan, 1, MPI_DOUBLE, 0, comm);

		// cout << "test_888: " << *tmp_xuan << endl;
		double sum_of_elements = 0.d;
		if(comRank == 0){
			for(int tmp=0; tmp < nProcess; tmp++){
				sum_of_elements += tmp_xuan[tmp];
			}
			cout << "test_999:  " << sum_of_elements << endl;
      string fname_xuan = "./example/Pressure_Drop.txt";
		  ofstream fout_xuan;
		  fout_xuan.open(fname_xuan);
		  fout_xuan << std::fixed << std::setprecision(8) << sum_of_elements << endl;
		  fout_xuan.close();
		} 

		MPI_Barrier(comm);
		/*Visualize the result*/
		if (comRank == 0){
			cout << "Visualizing...\n";
			VisualizeVTK_ControlMesh(cpts, tmesh, l, fn);
		}
		MPI_Barrier(comm);
		//VisualizeVTK_PhysicalDomain(l, fn + "final_physics");
	}
	MatDestroy(&GK);
	VecDestroy(&GR);
	VecDestroy(&temp_solution);
	KSPDestroy(&ksp);

	MPI_Barrier(comm);
}

void NS_2Dsteady::VisualizeVTK_ControlMesh(const vector<Vertex2D>& spt, const vector<Element2D>& mesh, int step, string fn)
{
	stringstream ss;
	ss << step;

	//string fname = fn + "controlmesh_VelocityPressure_"+ss.str()+".vtk";
	/*
  string fname = fn + "controlmesh_VelocityPressure.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i].coor[0] << " " << spt[i].coor[1] << " " << 0.000000 << "\n";
		}
		fout << "\nCELLS " << mesh.size() << " " << 5 * mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "4 " << mesh[i].IEN[0] << " " << mesh[i].IEN[1] << " " << mesh[i].IEN[2] << " " << mesh[i].IEN[3]
				<< " " << '\n';
		}
		fout << "\nCELL_TYPES " << mesh.size() << '\n';
		for (i = 0; i<mesh.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << Vel.size() / 2 << "\nVECTORS VelocityField float\n";
		for (uint i = 0; i<Vel.size() / 2; i++)
		{
			fout << Vel[i * 2] << " " << Vel[i * 2 + 1] << " " << 0.000000 << "\n";
		}
		fout << "\nSCALARS VelocityMagnitude float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<Pre.size(); i++)
		{
			fout << sqrt(Vel[i * 2] * Vel[i * 2] + Vel[i * 2 + 1] * Vel[i * 2 + 1]) << "\n";
		}
		fout << "\nSCALARS Pressure float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<Pre.size(); i++)
		{
			fout << Pre[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
  */
	string fname1 = fn + "velocityfield.txt";
	ofstream fout1;
	fout1.open(fname1.c_str());
	if (fout1.is_open())
	{
		for (uint i = 0; i<Vel.size() / 2; i++)
		{
			fout1 << Vel[i * 2] << " " << Vel[i * 2 + 1] << " " << 0.000000 /* << Vel[i * 3 + 2] */ << "\n";
		}
		fout1.close();
	}
	else
	{
		cout << "Cannot open " << fname1 << "!\n";
	}

}

void NS_2Dsteady::VisualizeVTK_PhysicalDomain(int step, string fn)
{
	vector<array<double, 3>> spt_all; // sample points
	vector<double> sresult_all;
	vector<array<int, 4>> sele_all;
	double detJ;
	int num_bzmesh_ele = bzmesh_process.size();
	double spt_proc[num_bzmesh_ele * 4 * 3];
	double sresult_proc[num_bzmesh_ele * 4 * 3];
	int sele_proc[num_bzmesh_ele * 4];

	for (unsigned int e = 0; e<num_bzmesh_ele; e++)
	{
		int ns(2);
		vector<double> su(ns);
		for (int i = 0; i < ns; i++)
		{
			su[i] = double(i) / (double(ns) - 1.);
		}

		int loc(0);
		for (int a = 0; a<ns; a++)
		{
			for (int b = 0; b<ns; b++)
			{
				double pt1[2], dudx[2];
				double result[3];
				ResultCal_Bezier(su[b], su[a], bzmesh_process[e], pt1, result, dudx, detJ);
				spt_proc[12 * e + loc*3 + 0] = pt1[0];
				spt_proc[12 * e + loc*3 + 1] = pt1[1];
				spt_proc[12 * e + loc*3 + 2] = 0.000000; // pt1[2];
				sresult_proc[12 * e + loc * 3 + 0] = result[0];
				sresult_proc[12 * e + loc * 3 + 1] = result[1];
				sresult_proc[12 * e + loc * 3 + 2] = result[2];
				// sresult_proc[32 * e + loc * 4 + 3] = result[3];
				loc++;
			}
		}
		int nns[2] = { ns*ns, ns };
		for (int a = 0; a < ns - 1; a++)
		{
			for (int b = 0; b < ns - 1; b++)
			{
				// for (int c = 0; c < ns - 1; c++)
				// {
				sele_proc[4 * e + 0] = 4 * e + a * ns + b;
				sele_proc[4 * e + 1] = 4 * e + a * ns + b + 1;
				sele_proc[4 * e + 2] = 4 * e + (a + 1) * ns + (b + 1);
				sele_proc[4 * e + 3] = 4 * e + (a + 1) * ns + b;
				// sele_proc[8 * e + 4] = 8 * e + (a + 1)*nns[1] + b*ns + c;
				// sele_proc[8 * e + 5] = 8 * e + (a + 1)*nns[1] + b*ns + c + 1;
				// sele_proc[8 * e + 6] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
				// sele_proc[8 * e + 7] = 8 * e + (a + 1)*nns[1] + (b + 1)*ns + c;
				// }
			}
		}
	}
	

	double *spts= NULL;
	double *sresults=NULL;
	int *seles=NULL;
	int *displs_spts = NULL;
	int *displs_sresults = NULL;
	int *displs_seles = NULL;
	int *num_bzmesh_eles=NULL;
	int *recvcounts_spts = NULL;
	int *recvcounts_sresults = NULL;
	int *recvcounts_seles = NULL;

	if (comRank == 0)
	{		
		num_bzmesh_eles = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_spts = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_sresults = (int*)malloc(sizeof(int)*nProcess);
		recvcounts_seles = (int*)malloc(sizeof(int)*nProcess);
	}
	MPI_Gather(&num_bzmesh_ele, 1, MPI_INT, num_bzmesh_eles, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Barrier(comm);

	if (comRank == 0)
	{
		spts = (double*)malloc(sizeof(double) * 12 * n_bzmesh);
		sresults = (double*)malloc(sizeof(double) * 12 * n_bzmesh);
		seles = (int*)malloc(sizeof(int) * 4 * n_bzmesh);

		displs_spts = (int*)malloc(nProcess * sizeof(int));
		displs_sresults = (int*)malloc(nProcess * sizeof(int));
		displs_seles = (int*)malloc(nProcess * sizeof(int));
		displs_spts[0] = 0;
		displs_sresults[0] = 0;
		displs_seles[0] = 0;

		for (int i = 1; i<nProcess; i++) 
		{
			displs_spts[i] = displs_spts[i - 1] + num_bzmesh_eles[i - 1] * 12;
			displs_sresults[i] = displs_sresults[i - 1] + num_bzmesh_eles[i - 1] * 12;
			displs_seles[i] = displs_seles[i - 1] + num_bzmesh_eles[i - 1] * 4;
		}

		for (int i = 0; i < nProcess; i++)
		{
			recvcounts_spts[i] = num_bzmesh_eles[i] * 12;
			recvcounts_sresults[i] = num_bzmesh_eles[i] * 12;
			recvcounts_seles[i] = num_bzmesh_eles[i] * 4;
		}
	}	

	MPI_Gatherv(spt_proc, num_bzmesh_ele * 4 * 3, MPI_DOUBLE, spts, recvcounts_spts, displs_spts, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sresult_proc, num_bzmesh_ele * 4 * 3, MPI_DOUBLE, sresults, recvcounts_sresults, displs_sresults ,MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	MPI_Gatherv(sele_proc, num_bzmesh_ele * 4, MPI_INT, seles, recvcounts_seles, displs_seles, MPI_INT, 0, PETSC_COMM_WORLD);
	
	
	if (comRank == 0)
	{
		for (int i = 0; i < n_bzmesh; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				array<double, 3> pt = { spts[i * 12 + j * 3 + 0], spts[i * 12 + j * 3 + 1], spts[i * 12 + j * 3 + 2] };
				spt_all.push_back(pt);
				sresult_all.push_back(sresults[i * 12 + j * 3 + 0]);
				sresult_all.push_back(sresults[i * 12 + j * 3 + 1]);
				sresult_all.push_back(sresults[i * 12 + j * 3 + 2]);
				// sresult_all.push_back(sresults[i * 32 + j * 4 + 3]);
			}
		}
		int sum_ele = 0;
		int pstart = 0;
		for (int i = 0; i < nProcess; i++)
		{
			for (int e = 0; e < num_bzmesh_eles[i]; e++)
			{
				array<int, 4> el;
				el[0] = pstart + seles[4 * sum_ele + 0];
				el[1] = pstart + seles[4 * sum_ele + 1];
				el[2] = pstart + seles[4 * sum_ele + 2];
				el[3] = pstart + seles[4 * sum_ele + 3];
				/*
				el[4] = pstart + seles[8 * sum_ele + 4];
				el[5] = pstart + seles[8 * sum_ele + 5];
				el[6] = pstart + seles[8 * sum_ele + 6];
				el[7] = pstart + seles[8 * sum_ele + 7];
				*/
				sele_all.push_back(el);
				sum_ele++;
			}
			pstart = pstart + num_bzmesh_eles[i] * 4;
		}
		cout << "Visualizing in Physical Domain...\n";
		WriteVTK(spt_all, sresult_all, sele_all, step, fn);
	}

}

void NS_2Dsteady::WriteVTK(const vector<array<double, 3>> spt, const vector<double> sdisp, const vector<array<int, 4>> sele, int step, string fn)
{
	string fname = fn + "_VelocityPressure.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
				<< " " /* << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] */ << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "9\n";
		}
		fout << "POINT_DATA " << sdisp.size() / 3 << "\nVECTORS VelocityField float\n";
		for (uint i = 0; i<sdisp.size() / 3; i++)
		{
			fout << sdisp[i * 3] << " " << sdisp[i * 3 + 1] << " " << 0.000000 /* << sdisp[i * 3 + 2] */<< "\n";
		}
		fout << "\nSCALARS VelocityMagnitude float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<sdisp.size() / 3; i++)
		{
			fout << sqrt(sdisp[i * 3] * sdisp[i * 3] + sdisp[i * 3 + 1] * sdisp[i * 3 + 1] /* + sdisp[i * 4 + 2] * sdisp[i * 4 + 2] */) << "\n";
		}
		fout << "\nSCALARS Pressure float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<sdisp.size() / 3; i++)
		{
			fout << sdisp[3 * i + 2] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void NS_2Dsteady::ResultCal_Bezier(double u, double v, const Element2D& bzel, double pt[2], double result[3], double dudx[3], double& detJ)
{
	double dUdx[2][2];
	vector<double> Nx(bzel.IEN.size());
	vector<array<double, 2>> dNdx(bzel.IEN.size());
	vector<array<array<double, 2>, 2>> dN2dx2;
	bzel.Para2Phys(u, v, pt);
	BasisFunction(u, v, bzel.IEN.size(), bzel.pts, bzel.cmat, Nx, dNdx, dN2dx2, dUdx, detJ);
	result[0] = 0.; result[1] = 0.; result[2] = 0.; // result[3] = 0.;
	for (uint i = 0; i < bzel.IEN.size(); i++)	{
		result[0] += Nx[i] * (Vel[dim* bzel.IEN[i] + 0]);
		result[1] += Nx[i] * (Vel[dim* bzel.IEN[i] + 1]);
		// result[2] += Nx[i] * (Vel[dim* bzel.IEN[i] + 2]);
		result[2] += Nx[i] * (Pre[bzel.IEN[i]]);
	}
}




