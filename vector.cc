// #pragma GCC optimize("unroll-loops", "omit-frame-pointer", "inline")
// //Optimization flags
#pragma GCC option("arch=native", "tune=native", "no-zero-upper") // Enable AVX
#pragma GCC target("avx")                                         // Enable
// AVX
#include <chrono>
#include <cstring>
#include <iostream>
#include <tuple>

#include <emmintrin.h>
#include <immintrin.h>

using namespace std;
using namespace std::chrono;

const int test_times = 100000;

const int world_size = sizeof(size_t);

void vector_copy(float *dest, float *src, int vector_size) {
  if (vector_size * sizeof(float) <
      world_size) { // byte copy assume word size 32/64
    *dest = *src;
    return;
  }
  int *d = (int *)dest;
  int *s = (int *)src;
  for (; vector_size >= world_size / sizeof(float) * 8;
       vector_size -= world_size / sizeof(float) * 8) {
    d[0] = s[0];
    d[1] = s[1];
    d[2] = s[2];
    d[3] = s[3];
    d[4] = s[4];
    d[5] = s[5];
    d[6] = s[6];
    d[7] = s[7];
    d += 8;
    s += 8;
  }
  for (; vector_size >= world_size / sizeof(float);
       vector_size -= world_size / sizeof(float)) {
    *d++ = *s++;
  }
  if (vector_size > 0)
    *(float *)d = *(float *)s;
  return;
}

void vector_mem_copy(float *dest, float *src, int vector_size) {
  memcpy(dest, src, vector_size * sizeof(float));
}
float vector_inner_product_naive(const float *a, const float *b,
                                 int vector_size) {

  float res = 0;
  for (int i = 0; i < vector_size; i++) {
    res += a[i] * b[i];
  }
  return res;
}
float vector_inner_product_simd(const float *a, const float *b,
                                int vector_size) {
  __m256 r0, r1, r2, r3;
  _mm256_zeroall();

  int i = 0;
  for (; i + 8 < vector_size; i += 8) {
    r1 = _mm256_loadu_ps(a + i);
    r2 = _mm256_loadu_ps(b + i);
    r3 = _mm256_mul_ps(r1, r2);
    r0 = _mm256_add_ps(r0, r3);
  }
  float res[8];
  _mm256_store_ps(res, r0);
  float ret = 0;
  for (int j = 0; j < 8; j++) {
    ret += res[i];
  }
  for (; i < vector_size; i++) {
    ret = a[i + 0] * b[i + 0];
  }
  return ret;
}

float vector_inner_product(const float *a, const float *b, int vector_size) {
  float res = 0;
  int i = 0;
  for (; i + 16 < vector_size; i += 16) {
    float res0 = a[i + 0] * b[i + 0];
    float res1 = a[i + 1] * b[i + 1];
    float res2 = a[i + 2] * b[i + 2];
    float res3 = a[i + 3] * b[i + 3];
    float res4 = a[i + 4] * b[i + 4];
    float res5 = a[i + 5] * b[i + 5];
    float res6 = a[i + 6] * b[i + 6];
    float res7 = a[i + 7] * b[i + 7];
    float res8 = a[i + 8] * b[i + 8];
    float res9 = a[i + 9] * b[i + 9];
    float res10 = a[i + 10] * b[i + 10];
    float res11 = a[i + 11] * b[i + 11];
    float res12 = a[i + 12] * b[i + 12];
    float res13 = a[i + 13] * b[i + 13];
    float res14 = a[i + 14] * b[i + 14];
    float res15 = a[i + 15] * b[i + 15];
    res0 = res0 + res1;
    res1 = res2 + res3;
    res2 = res4 + res5;
    res3 = res6 + res7;
    res4 = res8 + res9;
    res5 = res10 + res11;
    res6 = res12 + res13;
    res7 = res14 + res15;
    res0 = res0 + res1;
    res1 = res2 + res3;
    res2 = res4 + res5;
    res3 = res6 + res7;
    res0 = res0 + res1;
    res1 = res2 + res3;
    res += res0 + res1;
  }

  for (; i < vector_size; i++) {
    res += a[i] * b[i];
  }
  return res;
}

void vector_sum(float *c, float *a, float *b, int vector_size) {
  int i = 0;
  for (; i + 16 < vector_size; i += 16) {
    c[i + 0] = a[i + 0] + b[i + 0];
    c[i + 1] = a[i + 1] + b[i + 1];
    c[i + 2] = a[i + 2] + b[i + 2];
    c[i + 3] = a[i + 3] + b[i + 3];
    c[i + 4] = a[i + 4] + b[i + 4];
    c[i + 5] = a[i + 5] + b[i + 5];
    c[i + 6] = a[i + 6] + b[i + 6];
    c[i + 7] = a[i + 7] + b[i + 7];
    c[i + 8] = a[i + 8] + b[i + 8];
    c[i + 9] = a[i + 9] + b[i + 9];
    c[i + 10] = a[i + 10] + b[i + 10];
    c[i + 11] = a[i + 11] + b[i + 11];
    c[i + 12] = a[i + 12] + b[i + 12];
    c[i + 13] = a[i + 13] + b[i + 13];
    c[i + 14] = a[i + 14] + b[i + 14];
    c[i + 15] = a[i + 15] + b[i + 15];
  }
  for (; i < vector_size; i++) {
    c[i] = a[i] + b[i];
  }
}

void init_vector(float *v, int vector_size) {
  for (int i = 0; i < vector_size; i++) {
    v[i] = (float)i * 1.1234567f;
  }
}

void test_copy(int vector_size) {
  float *a = new float[vector_size];
  float *b = new float[vector_size];
  init_vector(b, vector_size);

  auto t0 = steady_clock::now();
  for (int i = 0; i < test_times; i++) {
    vector_copy(a, b, vector_size);
  }
  int nano_t = duration_cast<nanoseconds>(steady_clock::now() - t0).count();

  double flops = test_times * vector_size / (double)nano_t;
  cout << vector_size << " " << flops << endl;
}

void test_mem_copy(int vector_size) {
  float *a = new float[vector_size];
  float *b = new float[vector_size];
  init_vector(b, vector_size);

  auto t0 = steady_clock::now();
  for (int i = 0; i < test_times; i++) {
    vector_mem_copy(a, b, vector_size);
  }
  int nano_t = duration_cast<nanoseconds>(steady_clock::now() - t0).count();

  double flops = test_times * vector_size / (double)nano_t;
  cout << vector_size << " " << flops << endl;
}

void test_sum(int vector_size) {
  float *a = new float[vector_size];
  float *b = new float[vector_size];
  float *c = new float[vector_size];
  init_vector(a, vector_size);
  init_vector(b, vector_size);

  auto t0 = steady_clock::now();
  for (int i = 0; i < test_times; i++) {
    vector_sum(c, a, b, vector_size);
  }
  int nano_t = duration_cast<nanoseconds>(steady_clock::now() - t0).count();

  double flops = test_times * vector_size / (double)nano_t;
  cout << vector_size << " " << flops << endl;
}

static float what = 0;

void test_inner(int vector_size) {
  float *a = new float[vector_size];
  float *b = new float[vector_size];
  init_vector(a, vector_size);
  init_vector(b, vector_size);

  auto t0 = steady_clock::now();
  for (int i = 0; i < test_times; i++) {
    what = vector_inner_product_simd(a, b, vector_size);
  }
  int nano_t = duration_cast<nanoseconds>(steady_clock::now() - t0).count();

  double flops = test_times * vector_size / (double)nano_t;
  cout << vector_size << " " << flops * 2 << endl;
}

void test_all_length(void (*f)(int), const char *test_name) {

  cout << test_name << endl;
  for (int i = 2; i < 1026; i += 1) {
    f(i);
  }
  cout << "Done" << endl;
}

int main(int argc, const char *argv[]) {
  test_all_length(test_copy, "test copy");
  test_all_length(test_inner, "test inner");
  test_all_length(test_sum, "test sum");
  cout << what << endl;
  return 0;
}
