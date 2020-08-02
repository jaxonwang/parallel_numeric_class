#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <thread>
#include <unordered_map>
#include <vector>

#include <sys/types.h>
#include <unistd.h>

#include "datagen.h"

using namespace std;

static int world_size;
static int blk_size;

enum Tag {
  expected_solution,
  b_vector,
  matrix_block,
  L_matrix,
};

int get_rank_by_block_id(const int i, const int j) {
  // block-cyclic column partition
  return j % world_size;
}

void gen_and_distribute(const int len, const int blk_size) {
  assert(len % blk_size == 0);
  auto full_matrix = gen_matrix(len);
  auto solution = gen_solution(len);
  auto b = gen_b(full_matrix, solution);
  for (int i = 0; i < world_size; i++) {
    MPI_Send(solution.data(), solution.size(), MPI_DOUBLE, i,
             Tag::expected_solution, MPI_COMM_WORLD);
    MPI_Send(b.data(), b.size(), MPI_DOUBLE, i, Tag::b_vector, MPI_COMM_WORLD);
  }

  int block_divided_len = len / blk_size;
  for (int i = 0; i < block_divided_len; i++) {
    for (int j = 0; j <= i; j++) {
      vector<double> block_ij(blk_size * blk_size, 0.0);
      // copy a block
      for (int k = 0; k < blk_size; k++) {
        for (int l = 0; l < blk_size; l++) {
          block_ij[k * blk_size + l] =
              full_matrix[(i * blk_size + k) * len + j * blk_size + l];
        }
      }
      int target = get_rank_by_block_id(i, j);
      MPI_Send(block_ij.data(), block_ij.size(), MPI_DOUBLE, target,
               Tag::matrix_block, MPI_COMM_WORLD);
#ifdef DEBUG
      cout << "Sending block" << i << "," << j << " to " << target << endl;
#endif
    }
  }
}

namespace std {
template <> struct hash<pair<int, int>> {
  std::size_t operator()(const pair<int, int> &p) const noexcept {
    std::size_t h1 = std::hash<int>{}(p.first);
    std::size_t h2 = std::hash<int>{}(p.second);
    return h1 ^ (h2 << 1); // or use boost::hash_combine
  }
};
} // namespace std

typedef unordered_map<pair<int, int>, vector<double>> blocks_t;

vector<double> &get_blocks(blocks_t &blocks, const int i, const int j) {
  auto key = pair<int, int>{i, j};
  return blocks[key];
}

void receive_data(int rank, const int len, const int blk_size,
                  vector<double> &expected_solution, vector<double> &b_vector,
                  blocks_t &blocks) {
  expected_solution.resize(len);
  b_vector.resize(len);

  MPI_Recv(expected_solution.data(), len, MPI_DOUBLE, 0, Tag::expected_solution,
           MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  MPI_Recv(b_vector.data(), len, MPI_DOUBLE, 0, Tag::b_vector, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);

  int block_divided_len = len / blk_size;
  for (int i = 0; i < block_divided_len; i++) {
    for (int j = 0; j <= i; j++) {
      if (get_rank_by_block_id(i, j) != rank)
        continue;
      vector<double> block_ij(blk_size * blk_size, 0.0);
      // copy a block
      MPI_Recv(block_ij.data(), block_ij.size(), MPI_DOUBLE, 0,
               Tag::matrix_block, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // cout << rank << " received block " << i << "," << j << endl;
      blocks[pair<int, int>{i, j}] = move(block_ij);
    }
  }
}

void ldlt(const int rank, const int len, const int blk_size, blocks_t &blocks) {

  int block_divided_len = len / blk_size;
  for (int j = 0; j < block_divided_len; j++) {

    if (get_rank_by_block_id(0, j) == rank) {
      auto &Ajj = get_blocks(blocks, j, j);
      chol_block(Ajj, blk_size); // Ajj = chol(Ajj)
      auto Ajj_inverse = trangular_inverse(Ajj, blk_size);
      for (int i = j + 1; i < block_divided_len; i++) {
        auto &Aij = get_blocks(blocks, i, j);
        auto newAij =
            A_mul_Bt(Aij, Ajj_inverse, blk_size); // Aij = Aij * Ajj_-t
        Aij = move(newAij);
        for (int k = j + 1; k <= i; k++) {
          int target = get_rank_by_block_id(i, k);
          auto &Lkj = get_blocks(blocks, k, j);
          if (target == rank) { // block is local
            A_sub_B(get_blocks(blocks, i, k), A_mul_Bt(Aij, Lkj, blk_size),
                    blk_size); // Aik = Aik - Aij * Akj_t
          } else {
            MPI_Send(Aij.data(), Aij.size(), MPI_DOUBLE, target, Tag::L_matrix,
                     MPI_COMM_WORLD);
            MPI_Send(Lkj.data(), Lkj.size(), MPI_DOUBLE, target, Tag::L_matrix,
                     MPI_COMM_WORLD);
          }
        }
      }
    } else {
      for (int i = j + 1; i < block_divided_len; i++) {
        for (int k = j + 1; k <= i; k++) {
          int target = get_rank_by_block_id(i, k);
          if (target != rank)
            continue;
          int source = get_rank_by_block_id(i, j);
          vector<double> Aij(blk_size * blk_size, 0);
          vector<double> Lkj(blk_size * blk_size, 0);
          MPI_Recv(Aij.data(), Aij.size(), MPI_DOUBLE, source, Tag::L_matrix,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(Lkj.data(), Lkj.size(), MPI_DOUBLE, source, Tag::L_matrix,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          A_sub_B(get_blocks(blocks, i, k), A_mul_Bt(Aij, Lkj, blk_size),
                  blk_size); // Aik = Aik - Aij * Akj_t
        }
      }
    }
  }
}

void run(const int rank, const int len, const int blk_size) {
  unique_ptr<thread> t;
  if (rank == 0) {
    // start a new thread to distribute data
    t.reset(new thread{[=]() { gen_and_distribute(len, blk_size); }});
  }

  vector<double> expected_solution;
  vector<double> b_vector;
  unordered_map<pair<int, int>, vector<double>> blocks;
  receive_data(rank, len, blk_size, expected_solution, b_vector, blocks);

  ldlt(rank, len, blk_size, blocks);
  // cout << rank << "'s solution"
  //      << " ";
  // for (int i = 0; i < len; i++) {
  //   cout << expected_solution[i] << " ";
  // }
  // cout << endl;
  // cout << rank << "'s bvector"
  //      << " ";
  // for (int i = 0; i < len; i++) {
  //   cout << b_vector[i] << " ";
  // }
  // cout << endl;

  // print blocks one by one
  int block_divided_len = len / blk_size;
  for (int i = 0; i < block_divided_len; i++) {
    for (int j = 0; j <= i; j++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (get_rank_by_block_id(i, j) == rank) {
        auto &block = get_blocks(blocks, i, j);
        cout << rank << "'s block" << i << "," << j << endl;
        for (int k = 0; k < blk_size; k++) {
          for (int l = 0; l < blk_size; l++) {
            cout << block[k * blk_size + l] << " ";
          }
          cout << endl;
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  if (rank == 0) {
    t->join();
  }
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cout << "Usage: ldlt_mpi square_matrix_len block_size" << endl;
    return 1;
  }
  const int len = stoi(argv[1]);
  blk_size = stoi(argv[2]);
  if (len % blk_size != 0) {
    cout << "square_matrix_len should be divisible by blk_size" << endl;
    return 1;
  }
  int provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  if (provided < MPI_THREAD_MULTIPLE) {
    printf("ERROR: The MPI library does not have full thread support\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int my_world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_world_rank);

#ifdef DEBUG
  int ifl = 0;
  if (my_world_rank == 0) {
    cout << "PID " << getpid() << " is ready for attach." << endl;
    while (true) {
      if (!ifl) {
        continue;
      }
      cout << "PID " << getpid() << " is attached." << endl;
      break;
    }
  }
#endif
  run(my_world_rank, len, blk_size);

  MPI_Finalize();

  return 0;
}
