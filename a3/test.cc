#include <iostream>
#include <vector>

#include "datagen.h"

int main(int argc, char *argv[]) {

  int len = 20;
  auto matrix = gen_matrix(len);
  auto matrix_i = trangular_inverse(matrix, len);
  auto matrix_ii = trangular_inverse(matrix_i, len);
  print_matrix(matrix, len);
  print_matrix(matrix_i, len);
  print_matrix(matrix_ii, len);
  print_matrix(A_mul_B(matrix_ii, matrix_i, len), len);
  transpose(matrix_ii, len);
  print_matrix(A_mul_Bt(matrix_ii, matrix_i, len),len);

  return 0;
}
