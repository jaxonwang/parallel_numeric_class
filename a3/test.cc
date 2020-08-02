#include <iostream>
#include <vector>

#include "datagen.h"

bool equal(const vector<double> &A, const vector<double> &B, const int len){
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            if(A[i*len+j] - B[i*len+j] < 1e-15)
                continue;
            if(A[i*len+j] - B[i*len+j] < -1e-15)
                continue;
            return false;
        }
    }
    return true;
}

void assert(bool b){
    if(!b)
        cout << "assert failed" << endl;
    else{
        cout << "assert success" << endl;
    }
}

vector<double> gen_unit(const int len){
    vector<double> u(len * len, 0);
    for (int i = 0; i < len; i++) {
        u[i*len+i] = 1;
    }
    return u;
}

int main(int argc, char *argv[]) {

  int len = 20;
  auto unit_matrix = gen_unit(len);

  auto matrix = gen_matrix(len);
  auto matrix_i = trangular_inverse(matrix, len);
  auto matrix_ii = trangular_inverse(matrix_i, len);
  auto ret = A_mul_B(matrix_ii, matrix_i, len);
  assert(equal(ret, unit_matrix, len));
  transpose(matrix_ii, len);
  auto ret1 = A_mul_Bt(matrix_ii, matrix_i, len);
  assert(equal(ret1, unit_matrix, len));

  auto matrix_t = matrix;
  transpose(matrix_t, len);
  assert(!equal(matrix_t, matrix,len));
  transpose(matrix_t, len);
  assert(equal(matrix_t, matrix, len));

  return 0;
}
