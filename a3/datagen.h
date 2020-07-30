#pragma once

#include <vector>

using namespace std;

vector<double> symMatVec(vector<double> &a, vector<double> &x);
vector<double> gen_solution(const int len);
vector<double> gen_matrix(const int len);
vector<double> gen_b(vector<double> &a, vector<double> &x);
