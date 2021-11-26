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
using namespace std;

//problem setting
class UserSetting
{
public:
	UserSetting();
	void SetVariables(string fn_par, vector<double>& var);
	void SetInitialCondition(int ndof, vector<double>& Vel0, vector<double>& Pre0, vector<Vertex2D>& pts,  const vector<array<double, 2>> velocity_node, const int dim);
	void ReadMesh(string fn, vector<Vertex2D>& pts, vector<Element2D>& mesh);
	void ReadVelocityField(string fn, int npts, vector<array<double, 2>>& velocity);
	void AssignProcessor(string fn, int &n_bzmesh, vector<vector<int>> &ele_process);
};
#endif