#pragma once

#include <vector>

using namespace std;

void print_matrix(const vector<double> &A, const int len);
vector<double> symMatVec(vector<double> &a, vector<double> &x);
vector<double> trangular_inverse(const vector<double> &A, const int len);
void chol_block(vector<double> &a, const int len);
vector<double> A_mul_Bt(const vector<double> &A, const vector<double> &B, const int len);
vector<double> A_mul_B(const vector<double> &A, const vector<double> &B, const int len);
void transpose(vector<double>&A, const int len);
void A_sub_B(vector<double> &A, const vector<double> &B, const int len);
vector<double> gen_solution(const int len);
vector<double> gen_matrix(const int len);
vector<double> gen_b(vector<double> &a, vector<double> &x);
