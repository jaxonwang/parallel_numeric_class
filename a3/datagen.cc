#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void print_matrix(const vector<double> &A, const int len) {
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            cout << A[i * len + j] << " ";
        }
        cout << endl;
    }
}

vector<double> symMatVec(vector<double> &a, vector<double> &x) {
    vector<double> y(x.size(), 0);
    int i, j;
    int n = x.size();
    for (i = 0; i < n; i++) {
        double t = 0.0;
        for (j = 0; j <= i; j++) t += a[i * n + j] * x[j];

        for (j = i + 1; j < n; j++) t += a[j * n + i] * x[j];

        y[i] = t;
    }
    return y;
}

vector<double> trangular_inverse(const vector<double> &A, const int len) {
    vector<double> tmp(A.size(), 0);
    for (int j = 0; j < len; j++) {
        tmp[j * len + j] = 1.0 / A[j * len + j];
        for (int i = j + 1; i < len; i++) {
            double s = 0.0;
            for (int k = j; k < i; k++) {
                s += A[i * len + k] * tmp[k * len + j];
            }
            tmp[i * len + j] = -s / A[i * len + i];
        }
    }
    return tmp;
}

void chol_block(vector<double> &a, const int len) {
    int n = len;
    for (int i = 0; i < n; i++) {
        a[i * n + i] = sqrt(a[i * n + i]);
        double invp = 1.0 / a[i * n + i];
        for (int j = i + 1; j < n; j++) {
            a[j * n + i] *= invp;
            for (int k = i + 1; k <= j; k++)
                a[j * n + k] -= a[j * n + i] * a[k * n + i];
        }
    }
}

vector<double> A_mul_Bt(const vector<double> &A, const vector<double> &B,
                        const int len) {
    vector<double> C(B.size(), 0);
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            double tmp = 0;
            for (int k = 0; k < len; k++) {
                tmp += A[i * len + k] * B[j * len + k];
            }
            C[i * len + j] = tmp;
        }
    }
    return C;
}

void transpose(vector<double> &A, const int len) {
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < i; j++) {
            swap(A[i * len + j], A[j * len + i]);
        }
    }
}

vector<double> A_mul_B(const vector<double> &A, const vector<double> &B,
                       const int len) {
    vector<double> C(B.size(), 0);
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            double tmp = 0;
            for (int k = 0; k < len; k++) {
                tmp += A[i * len + k] * B[k * len + j];
            }
            C[i * len + j] = tmp;
        }
    }
    return C;
}

void A_sub_B(vector<double> &A, const vector<double> &B, const int len) {
    vector<double> C(B.size(), 0);
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            A[i * len + j] -= B[i * len + j];
        }
    }
}

vector<double> gen_solution(const int len) { return vector<double>(len, 1.0); }

vector<double> gen_matrix(const int len) {
    vector<double> a(len * len, 0);

    int i, j;
    int n = len;

    /* fill lower triangular elements */
    for (i = 0; i < n; i++)
        for (j = 0; j < i; j++) a[i * n + j] = rand() / (RAND_MAX + 1.0);

    /* fill diagonal elements */
    for (i = 0; i < n; i++) {
        double s = 0.0;
        for (j = 0; j < i; j++) s += a[i * n + j];

        for (j = i + 1; j < n; j++) s += a[j * n + i]; /* upper triangular */

        a[i * n + i] = s + 1.0; /* diagonal dominant */
    }
    return a;
}

vector<double> gen_b(vector<double> &a, vector<double> &x) {
    return symMatVec(a, x);
}

void gen_data(const int len, const int blk_size, const string &dir) {
    ofstream f_m(dir + "/full_matrix.data", ios_base::trunc | ios_base::out);
    ofstream f_b(dir + "/b.data", ios_base::trunc | ios_base::out);
    ofstream f_s(dir + "/solution.data", ios_base::trunc | ios_base::out);
    auto full_matrix = gen_matrix(len);
    auto solution = gen_solution(len);
    auto b = gen_b(full_matrix, solution);

    auto vector_writer = [](ofstream &fs, vector<double> v) {
        if (!v.size()) return;
        size_t i = 0;
        for (; i < v.size() - 1; i++) {
            fs << v[i] << endl;
        }
        fs << v[i];
    };

    auto matrix_writer = [](ofstream &fs, vector<double> matrix,
                            const int len) {
        for (int i = 0; i < len; i++) {
            int j = 0;
            for (; j < len - 1; j++) {
                fs << matrix[i * len + j] << " ";
            }
            fs << matrix[j] << endl;
        }
    };

    vector_writer(f_b, b);
    vector_writer(f_s, solution);

    matrix_writer(f_m, full_matrix, len);

    int block_divided_len = len / blk_size;

    for (int i = 0; i < block_divided_len; i++) {
        for (int j = 0; j <= i; j++) {
            vector<double> block_ij(blk_size * blk_size, 0.0);
            ofstream block_fs{dir + "/A" + to_string(i) + to_string(j),
                              ios_base::out | ios_base::trunc};
            // copy a block
            for (int k = 0; k < blk_size; k++) {
                for (int l = 0; l < blk_size; l++) {
                    block_ij[k * blk_size + l] =
                        full_matrix[(i * blk_size + k) * len + j * blk_size +
                                    l];
                }
            }
            matrix_writer(block_fs, block_ij, blk_size);
        }
    }
}

// int main(int argc, char *argv[]) {
//   if (argc != 4) {
//     cout << "Usage: a.out square_matrix_len block_size dst_dir" << endl;
//     return 1;
//   }
//
//   const int len = stoi(argv[1]);
//   const int blk_size = stoi(argv[2]);
//   if (len % blk_size != 0) {
//     cout << "square_matrix_len should be divisible by blk_size" << endl;
//     return 1;
//   }
//   const string dst_dir{argv[3]};
//
//   gen_data(len, blk_size, dst_dir);
//
//   return 0;
// }
