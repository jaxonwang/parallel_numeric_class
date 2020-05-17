#pragma GCC optimize("unroll-loops", "omit-frame-pointer", "inline")
// //Optimization flags
#pragma GCC option("arch=native", "tune=native", "no-zero-upper") // Enable AVX
#pragma GCC target("avx")                                         // Enable
// AVX
#include <chrono>
#include <iostream>
#include <tuple>

using namespace std;
using namespace std::chrono;

const int test_times = 10000;

template <int vector_size> void vector_copy(float *dest, float *src) {
  for (int i = 0; i < vector_size; i++) {
    dest[i] = src[i];
  }
}

template <int vector_size>
void vector_inner_product(float *c, float *a, float *b) {
  for (int i = 0; i < vector_size; i++) {
    c[i] = a[i] * b[i];
  }
}

template <int vector_size> void vector_sum(float *c, float *a, float *b) {
  for (int i = 0; i < vector_size; i++) {
    c[i] = a[i] + b[i];
  }
}

template <int vector_size> void init_vector(float *v) {
  for (int i = 0; i < vector_size; i++) {
    v[i] = (float)i * 1.1234567f;
  }
}

template <int vector_size> void test_copy() {
  float *a = new float[vector_size];
  float *b = new float[vector_size];
  init_vector<vector_size>(b);

  auto t0 = steady_clock::now();
  for (int i = 0; i < test_times; i++) {
    vector_copy<vector_size>(a, b);
  }
  int nano_t = duration_cast<nanoseconds>(steady_clock::now() - t0).count();

  double flops = test_times * vector_size / (double)nano_t;
  cout << "Test vector copy with size: " << vector_size << endl;
  cout << "GFLOPS: " << flops << endl;
}

template <int vector_size> void test_compute() {
  float *a = new float[vector_size];
  float *b = new float[vector_size];
  float *c = new float[vector_size];
  init_vector<vector_size>(a);
  init_vector<vector_size>(b);

  auto t0 = steady_clock::now();
  for (int i = 0; i < test_times; i++) {
    vector_inner_product<vector_size>(c, a, b);
  }
  int nano_t = duration_cast<nanoseconds>(steady_clock::now() - t0).count();

  double flops = test_times * vector_size / (double)nano_t;
  cout << "Test vector inner product with size: " << vector_size << endl;
  cout << "GFLOPS: " << flops << endl;
}

template <int N = 2> constexpr auto gen_len_list() {
  auto t = make_tuple(N);
  if constexpr (N == 32) {
    return t;
  } else {
    return tuple_cat(t, gen_len_list<N + 1>());
  }
}
constexpr auto lengh_list = gen_len_list();

constexpr auto lengh_list = make_tuple(
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
    40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
    59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
    78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96,
    97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112,
    113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
    128);

template <int N = 0> void gen_test_copy() {
  test_copy<get<N>(lengh_list)>();
  if constexpr (N + 1 < tuple_size<decltype(lengh_list)>::value)
    gen_test_copy<N + 1>();
}

template <int N = 0> void gen_test_compute() {
  test_compute<get<N>(lengh_list)>();
  if constexpr (N + 1 < tuple_size<decltype(lengh_list)>::value)
    gen_test_compute<N + 1>();
}

int main(int argc, const char *argv[]) {
  gen_test_copy();
  gen_test_compute();
  return 0;
}
