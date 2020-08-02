#include <cmath>
#include <iostream>

#include "datagen.h"

vector<double> solveSym(vector<double> &a, vector<double> &b) {
  int i, j, k;
  int n = b.size();

  vector<double> x(b.size(), 0);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      cout << a[i * n+ j] << " ";
    }
    cout << endl;
  }

  /* LDLT decomposition: A = L * D * L^t */
  for (i = 0; i < n; i++) {
    double invp = 1.0 / a[i * n + i];

    for (j = i + 1; j < n; j++) {
      double aji = a[j * n + i];
      a[j * n + i] *= invp;

      for (k = i + 1; k <= j; k++)
        a[j * n + k] -= aji * a[k * n + i];
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      cout << a[i * n+ j] << " ";
    }
    cout << endl;
  }

  /* forward solve L y = b: but y is stored in x
     can be merged to the previous loop */
  for (i = 0; i < n; i++) {
    double t = b[i];

    for (j = 0; j < i; j++)
      t -= a[i * n + j] * x[j];

    x[i] = t;
  }

  /* backward solve D L^t x = y */
  for (i = n - 1; i >= 0; i--) {
    double t = x[i] / a[i * n + i];

    for (j = i + 1; j < n; j++)
      t -= a[j * n + i] * x[j];

    x[i] = t;
  }

  return x;
}

int main(int argc, char *argv[]) {

  int n, i;
  if (argc != 2) {
    cout << "Usage: ldlt_mpi square_matrix_len" << endl;
    return 1;
  }
  n = stoi(argv[1]);

  auto a = gen_matrix(n);
  auto solution = gen_solution(n);
  auto b = gen_b(a, solution);

  auto x = solveSym(a, b);
  double e = 0;
  for (i = 0; i < n; i++)
    e += (x[i] - solution[i]) * (x[i] - solution[i]);
  e = sqrt(e);
  cout << "error norm = " << e << endl;
  cout << "--- good if error is around n * 1e-16 or less" << endl;

  return 0;
}
