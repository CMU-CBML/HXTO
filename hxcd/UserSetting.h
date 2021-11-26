#ifndef USERSETTING_H
#define USERSETTING_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "BasicDataStructure.h"
#include <cmath>
#ifdef _WIN32
#include <io.h>
#include <direct.h>
#endif
#ifdef __linux__ 
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif
using namespace std;

//problem setting
class UserSetting
{
public:
	UserSetting();
	void SetVariables(string fn_par, vector<double>& var, string &path_out);
	void SetInitialCondition(vector<Vertex2D>& pts, vector<double>& var, vector<double>& CA0, vector<double>& NX0, vector<double>& NB0);
	void ReadMesh(string fn, vector<Vertex2D>& pts, vector<Element2D>& mesh);
	void ReadVelocityField(string fn, int npts, vector<array<double, 2>>& velocity);
	void AssignProcessor(string fn, int &n_bzmesh, vector<vector<int>> &ele_process);
};
#endif