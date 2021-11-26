#include "BasicDataStructure.h"

Vertex2D::Vertex2D()
{
	coor[0] = 0.;	coor[1] = 0.;	// coor[2] = 0.;
	label = 0;
}

Element2D::Element2D(int p)
{
	degree = p;
	order = p + 1;
	nbf = order*order;
	IEN.resize(4);
	pts.resize(4);
	for (int i = 0; i < 4; i++)
	{
		IEN[i] = 0;
		pts[i][0] = 0.; pts[i][1] = 0.; // pts[i][2] = 0.;
	}
}

void Element2D::BezierPolyn(double u, vector<double>& Nu, vector<double>& dNdu) const
{
	if (degree == 3)
	{
		double Nu0[4] = { (1. - u)*(1. - u)*(1. - u),3.*(1. - u)*(1. - u)*u,3.*(1. - u)*u*u,u*u*u };
		double dNdu0[4] = { -3.*(1. - u)*(1. - u),3. - 12.*u + 9.*u*u,3.*(2. - 3.*u)*u,3.*u*u };
		Nu.resize(order);
		dNdu.resize(order);
		for (int i = 0; i<order; i++)
		{
			Nu[i] = Nu0[i];
			dNdu[i] = dNdu0[i];
		}
	}
	else if (degree == 4)
	{
		double Nu0[5] = { (1. - u)*(1. - u)*(1. - u)*(1. - u),4.*(1. - u)*(1. - u)*(1. - u)*u,6.*(1. - u)*(1. - u)*u*u,4.*(1. - u)*u*u*u,u*u*u*u };
		double dNdu0[5] = { -4.*(1. - u)*(1. - u)*(1. - u),4.*(1. - u)*(1. - u)*(1. - 4.*u),12.*u*(1. - 3.*u + 2.*u*u),4.*(3. - 4.*u)*u*u,4.*u*u*u };
		Nu.resize(order);
		dNdu.resize(order);
		for (int i = 0; i<order; i++)
		{
			Nu[i] = Nu0[i];
			dNdu[i] = dNdu0[i];
		}
	}
}

void Element2D::Basis(double u, double v, vector<double>& Nt, vector<array<double, 2>>& dNdt) const
{
	vector<double> Nu, Nv, dNdu, dNdv;
	BezierPolyn(u, Nu, dNdu);
	BezierPolyn(v, Nv, dNdv);
	// BezierPolyn(w, Nw, dNdw);
	Nt.resize(nbf);
	dNdt.resize(nbf);
	int i, j, loc(0);
	for (j = 0; j<order; j++)
	{
		for (i = 0; i<order; i++)
		{
			// for (i = 0; i<order; i++)
			// {
				Nt[loc] = Nu[i] * Nv[j]; // * Nw[k];
				dNdt[loc][0] = dNdu[i] * Nv[j]; // * Nw[k];
				dNdt[loc][1] = Nu[i] * dNdv[j]; // * Nw[k];
				// dNdt[loc][2] = Nu[i] * Nv[j] * dNdw[k];
				loc++;
			// }
		}
	}
}

void Element2D::Para2Phys(double u, double v, double pt[2]) const
{
	vector<double> Nt;
	vector<array<double, 2>> dNdt;
	Basis(u, v, Nt, dNdt);
	pt[0] = 0.; pt[1] = 0.; // pt[2] = 0.;
	for (int i = 0; i<nbf; i++)
	{
		pt[0] += pts[i][0] * Nt[i];
		pt[1] += pts[i][1] * Nt[i];
		// pt[2] += pts[i][2] * Nt[i];
	}
}

void Raw2Vtk_hex(string fn)
{
	unsigned int npt, nel;
	vector<array<double, 2>> pts;
	vector<array<int, 4>> cnct;
	double tmp;
	string fn1(fn + ".raw");
	ifstream fin;
	fin.open(fn1);
	if (fin.is_open())
	{
		fin >> npt >> nel;
		pts.resize(npt);
		cnct.resize(nel);
		for (unsigned int i = 0; i < npt; i++)
		{
			fin >> pts[i][0] >> pts[i][1] >> tmp >> tmp;
			fin >> tmp;
		}
		for (unsigned int i = 0; i < nel; i++)
		{
			for (int j = 0; j < 4; j++) fin >> cnct[i][j];
			fin >> tmp;
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fn1 << "!\n";
	}
	string fn2(fn + ".vtk");
	ofstream fout;
	fout.open(fn2);
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		for (unsigned int i = 0; i<pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << 0.000000 /* << pts[i][2] */ << "\n";
		}
		fout << "\nCELLS " << cnct.size() << " " << 5 * cnct.size() << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			fout << "4 ";
			for (int j = 0; j<4; j++)
			{
				fout << cnct[i][j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << cnct.size() << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			fout << "9\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Cannot open " << fn2 << "!\n";
	}
}

void Rawn2Vtk_hex(string fn)
{
	unsigned int npt, nel;
	vector<array<double, 2>> pts;
	vector<array<int, 4>> cnct;
	double tmp;
	string fn1(fn + ".rawn");
	ifstream fin;
	fin.open(fn1);
	if (fin.is_open())
	{
		fin >> npt >> nel;
		pts.resize(npt);
		cnct.resize(nel);
		for (unsigned int i = 0; i < npt; i++)
		{
			fin >> pts[i][0] >> pts[i][1] >> tmp >> tmp >> tmp >> tmp >> tmp;
		}
		for (unsigned int i = 0; i < nel; i++)
		{
			for (int j = 0; j < 4; j++) fin >> cnct[i][j];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fn1 << "!\n";
	}
	string fn2(fn + ".vtk");
	ofstream fout;
	fout.open(fn2);
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		for (unsigned int i = 0; i<pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << 0.000000 /* pts[i][2] */ << "\n";
		}
		fout << "\nCELLS " << cnct.size() << " " << 5 * cnct.size() << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			fout << "4 ";
			for (int j = 0; j<4; j++)
			{
				fout << cnct[i][j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << cnct.size() << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			fout << "9\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Cannot open " << fn2 << "!\n";
	}
}

