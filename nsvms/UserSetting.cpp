#include "UserSetting.h"

UserSetting::UserSetting()
{
}

void UserSetting::SetVariables(string fn_par, vector<double>& var)
{
	var.resize(12);
	string fname(fn_par), stmp;
	stringstream ss;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 12; i++)
		{
			fin >> stmp >> var[i];
		}	
		fin.close();
	}	
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
	//double var0[2] = { 1.0,0.5 };	//nu, rou
	//var.clear();
	//var.resize(2);
	//var.assign(var0, var0 + 2);
}

void UserSetting::SetInitialCondition(int ndof, vector<double>& Vel0, vector<double>& Pre0, vector<Vertex2D>& pts, const vector<array<double, 2>> velocity_node, const int dim)
{
	double val_ini[2] = { 0.0, 0.0 }; //Vel0, Pre0
	Vel0.clear();
	Pre0.clear();
	Vel0.resize(ndof * dim, val_ini[0]);
	Pre0.resize(ndof, val_ini[1]);
}

void UserSetting::ReadMesh(string fn, vector<Vertex2D>& pts, vector<Element2D>& mesh)
{
	string fname(fn), stmp;
	int npts, neles, itmp;
	double tmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i < 4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		for (int i = 0; i < npts; i++)
		{
			fin >> pts[i].coor[0] >> pts[i].coor[1] >> tmp; // >> pts[i].coor[2]
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		mesh.resize(neles);
		for (int i = 0; i < neles; i++)
		{
			fin >> itmp >> mesh[i].IEN[0] >> mesh[i].IEN[1] >> mesh[i].IEN[2] >> mesh[i].IEN[3]; /* >>
				mesh[i].IEN[4] >> mesh[i].IEN[5] >> mesh[i].IEN[6] >> mesh[i].IEN[7]; */
			for (int j = 0; j < 4; j++)
			{
				mesh[i].pts[j][0] = pts[mesh[i].IEN[j]].coor[0];
				mesh[i].pts[j][1] = pts[mesh[i].IEN[j]].coor[1];
				// mesh[i].pts[j][2] = pts[mesh[i].IEN[j]].coor[2];
			}

		}
		for (int i = 0; i < neles + 5; i++) { getline(fin, stmp); } //skip lines

		string fname_xuan = "./example/Labels.txt";
		ifstream fin_xuan;
	  fin_xuan.open(fname_xuan);
		for (int i = 0; i < npts; i++)	{fin_xuan >> pts[i].label;}
	  fin_xuan.close();

		//for (int i = 0; i < npts; i++)	{ fin >> pts[i].label; }
		fin.close();
		//PetscPrintf(PETSC_COMM_WORLD, "Mesh Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}

}

void UserSetting::ReadVelocityField(string fn, int npts, vector<array<double, 2>>& velocity)
{
	string fname(fn);
	ifstream fin;
	fin.open(fname);
	velocity.resize(npts);
	double tmp;
	if (fin.is_open())
	{
		for (int i = 0; i < npts; i++)
		{
			fin >> velocity[i][0] >> velocity[i][1]; // >> tmp; // >> velocity[i][2];
		}
		fin.close();
    // cout << "xuan_March_01:   " << velocity[199][0] << endl;
		//PetscPrintf(PETSC_COMM_WORLD, "Velocity Field Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

void UserSetting::AssignProcessor(string fn, int &n_bzmesh, vector<vector<int>> &ele_process)
{
	int tmp;
	int i = 0;
	string fname(fn);
	ifstream fin;
	fin.open(fname, ios::in);
	if (fin.is_open())
	{
		while(!fin.eof() && fin.peek() != EOF)
		{
			fin >> tmp;
			ele_process[tmp].push_back(i);
			i++;
			fin.get();
		}
		n_bzmesh = i;
		//PetscPrintf(PETSC_COMM_WORLD, "Mesh partition finished!\n");
		fin.close();
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}
