#include "UserSetting.h"

UserSetting::UserSetting()
{
}

void UserSetting::SetVariables(string fn_par, vector<double>& var, string &path_out)
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
			if (i < 11)
			{
				ss << stmp << var[i] << "_";
			}				
			else 
			{
				ss << stmp << var[i];
			}
				
		}
		path_out = path_out + ss.str();
		
		fin.close();
		//PetscPrintf(PETSC_COMM_WORLD, "Parameter Loaded!\n");
		if (access(path_out.c_str(), 0) == -1)
		{
			PetscPrintf(PETSC_COMM_WORLD, "%s is not existing\n", path_out.c_str());
	
#ifdef _WIN32
			int flag = mkdir(path_out.c_str());
#endif
#ifdef __linux__ 
			int flag = mkdir(path_out.c_str(), 0777);
#endif
			if (flag == 0)
			{
				PetscPrintf(PETSC_COMM_WORLD, "%s make successfully\n", path_out.c_str());
			}
			else 
			{
				PetscPrintf(PETSC_COMM_WORLD, "%s make errorly\n", path_out.c_str());
			}
		}
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
	
	
	//double var0[7] = { 1.0,1.0,-0., //Dn0, v_plus, v_minus
	//					1.0,0.0,	// k+, k-
	//					0.1,0.0	};	//k'+,k'-
	//var.clear();
	//var.resize(7);
	//var.assign(var0, var0 + 7);
}

void UserSetting::SetInitialCondition(vector<Vertex2D>& pts, vector<double>& var, vector<double>& N0_ini, vector<double>& Nplus_ini, vector<double>& Nminus_ini)
{
	double val_ini[2] = { 0.0, 0.0 };//CA, N_plus, N_minus
	int npts = pts.size();
	N0_ini.clear();
	Nplus_ini.clear();
	Nminus_ini.clear();
	N0_ini.resize(npts, val_ini[0]);
	Nplus_ini.resize(npts, val_ini[1]);
	// Nminus_ini.resize(npts, val_ini[2]);

    // cout << npts << endl;
	for (int i = 0; i < npts; i++)
	{
		/// Default initial condition
		if (pts[i].label == 1)
		{
			N0_ini[i] = var[9];
			Nplus_ini[i] = var[10];
			// Nminus_ini[i] = var[11];
		}

		/// Initial condition for photoactivation experiment verification
		////movie2 photoactivation ICs
		//if (pts[i].coor[0] <17.5)
		//{
		//	CA0[i] = 1.0 * pow(abs(pts[i].coor[0] / 17.5), 1.3);
		//	NX0[i] = 2.0 * pow(abs(pts[i].coor[0] / 17.5), 1.3);
		//}
		//if (pts[i].coor[0]>=17.5 && pts[i].coor[0]<18.5)
		//{
		//	CA0[i] = 1.0;
		//	NX0[i] = 2.0;
		//}
		//if (pts[i].coor[0] > 18.5)
		//{
		//	CA0[i] = 1.0 - 0.85 / 15.5*(pts[i].coor[0] - 18.5);
		//	NX0[i] = 2.0 - 1.70 / 15.5*(pts[i].coor[0] - 18.5);
		//	//NX0[i] = 1.0 - 0.85 / 17*(pts[i].coor[0] - 17);
		//	//CA0[i] = 2.0 - 1.70 / 17*(pts[i].coor[0] - 17);
		//}

		////movie5 photoactivation ICs
		//if (pts[i].coor[0] <16.5)
		//{
		//	CA0[i] = 1.0 * pow(pts[i].coor[0] / 16.5, 2.0);
		//	NX0[i] = 2.0 * pow(pts[i].coor[0] / 16.5, 2.0);
		//}
		//if (pts[i].coor[0]>16.5 && pts[i].coor[0]<17.5)
		//{
		//	CA0[i] = 1.0;
		//	NX0[i] = 2.0;
		//}
		//if (pts[i].coor[0] > 17.5 && pts[i].coor[0]<25.5)
		//{
		//	CA0[i] = 1.0  * pow((31 - pts[i].coor[0]) / 13.5, 1.5);
		//	NX0[i] = 2.0  * pow((31 - pts[i].coor[0]) / 13.5, 1.5);
		//	//NX0[i] = 1.0 - 0.85 / 17*(pts[i].coor[0] - 17);
		//	//CA0[i] = 2.0 - 1.70 / 17*(pts[i].coor[0] - 17);
		//}
		//if (pts[i].coor[0] > 25.5 && pts[i].coor[0]<31.0 && pts[i].coor[1] <= -3.0)
		//{
		//	CA0[i] = 1.0  * pow((31 - pts[i].coor[0]) / 13.5, 1.5);
		//	NX0[i] = 2.0  * pow((31 - pts[i].coor[0]) / 13.5, 1.5);
		//	//NX0[i] = 1.0 - 0.85 / 17*(pts[i].coor[0] - 17);
		//	//CA0[i] = 2.0 - 1.70 / 17*(pts[i].coor[0] - 17);
		//}
		//if (pts[i].coor[0] > 25.5 && pts[i].coor[0]<32.5 && pts[i].coor[1]>-3.0)
		//{
		//	CA0[i] = 0.26*1.0* pow((32.5 - pts[i].coor[0]) / 7.0, 2.5);
		//	NX0[i] = 0.26*2.0* pow((32.5 - pts[i].coor[0]) / 7.0, 2.5);
		//	//NX0[i] = 1.0 - 0.85 / 17*(pts[i].coor[0] - 17);
		//	//CA0[i] = 2.0 - 1.70 / 17*(pts[i].coor[0] - 17);
		//}
	}
	//PetscPrintf(PETSC_COMM_WORLD, "User Set Initial Condition Ended!\n");
}

void UserSetting::ReadMesh(string fn, vector<Vertex2D>& pts, vector<Element2D>& mesh)//need vtk file with point label
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
			fin >> pts[i].coor[0] >> pts[i].coor[1] >> tmp; // pts[i].coor[2];
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
		for (int i = 0; i < neles + 5; i++) getline(fin, stmp);//skip lines
		string fname_xuan = "./example/Labels.txt";
		ifstream fin_xuan;
		fin_xuan.open(fname_xuan);
		for (int i = 0; i < npts; i++) {	fin_xuan >> pts[i].label; }
		fin_xuan.close();

		fin.close();
		//PetscPrintf(PETSC_COMM_WORLD, "Mesh Loaded!\n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Cannot open %s!\n", fname.c_str());
	}
}

void UserSetting::ReadVelocityField(string fn, int npts, vector<array<double,2>>& velocity)
{
	string fname(fn);
	ifstream fin;
	fin.open(fname);
	double tmp;
	velocity.resize(npts);
	if (fin.is_open())
	{
		for (int i = 0; i < npts; i++)
		{
			fin >> velocity[i][0] >> velocity[i][1] >> tmp; // >> velocity[i][2];
		}
		fin.close();
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
		while (!fin.eof() && fin.peek() != EOF)
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
