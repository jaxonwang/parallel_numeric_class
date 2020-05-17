#pragma GCC optimize("O3", "unroll-loops", "omit-frame-pointer", "inline")
// //Optimization flags
#pragma GCC option("arch=native", "tune=native", "no-zero-upper") // Enable AVX
#pragma GCC target("avx")                                         // Enable
// AVX
#include <chrono>
#include <iostream>
// #include <x86intrin.h> //AVX/SSEExtensions

using namespace std;
using namespace std::chrono;

int main(int argc, const char *argv[]) {

  const long int array_length = 512;
  const long int test_number = 100*1000;

  auto t0 = steady_clock::now();
  const float a = 0.000001234f;
  register float result0 = 0.00000024f;
  register float result1 = 0.00000023f;
  register float result2 = 0.00000022f;
  register float result3 = 0.00000021f;
  register float result4 = 0.00000020f;
  register float result5 = 0.00000019f;
  register float result6 = 0.00000018f;
  register float result7 = 0.00000017f;
  register float result8 = 0.00000016f;
  register float result9 = 0.00000015f;
  register float resulta = 0.00000014f;
  register float resultb = 0.00000013f;
  register float resultc = 0.00000012f;
  register float resultd = 0.00000011f;
  register float resulte = 0.00000010f;
  register float resultf = 0.00000009f;
  // register float result2 = 0;
  for (int i = 0; i < test_number; i++) {
    for (int i = 0; i < array_length; i++) {
        result0 += result0 * a;
        // result1 += result1 * a;
        // result2 += result2 * a;
        // result3 += result3 * a;
        // result4 += result4 * a;
        // result5 += result5 * a;
        // result6 += result6 * a;
        // result7 += result7 * a;
        // result8 += result8 * a;
        // result9 += result8 * a;
        // resulta += resulta * a;
        // resultb += resultb * a;
        // resultc += resultc * a;
        // resultd += resultd * a;
    }
  }

  auto t1 = steady_clock::now();
  // in nano second
  auto nano_t = duration_cast<nanoseconds>(t1 - t0).count();

  double flops = test_number * array_length / (double)nano_t;
  const int ops = 28;

  cout << "duration: " << nano_t << " ns." << endl;
  cout << "FLOPS: " << flops * ops << " GFLOPS" << endl;
  cout << "resulst " << result1 << " " << result2<< endl;
  cout << "resulst " << result3 << " " << result4<< endl;
  cout << "resulst " << result5 << " " << result6<< endl;
  cout << "resulst " << result7 << " " << result8<< endl;
  cout << "resulst " << result9 << " " << result0<< endl;
  cout << "resulst " << resulta << " " << resultb<< endl;
  cout << "resulst " << resultc << " " << resultd<< endl;

  return 0;
}
