#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>

#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <thread>
#include <unordered_map>
#include <vector>

#include "datagen.h"

using namespace std;

static int world_size;
static int blk_size;

enum Tag {
    expected_solution,
    b_vector,
    matrix_block,
    L_matrix,
    collect_result,
};

int get_rank_by_block_id(const int i, const int j) {
#ifdef LEFTLOOKING
    // block-cyclic row partition for right looking
    return i % world_size;
#else
    // block-cyclic column partition for left looking
    return j % world_size;
#endif
}

void gen_and_distribute(const int len, const int blk_size) {
    assert(len % blk_size == 0);
    auto full_matrix = gen_matrix(len);
    auto solution = gen_solution(len);
    auto b = gen_b(full_matrix, solution);
    for (int i = 0; i < world_size; i++) {
        MPI_Send(solution.data(), solution.size(), MPI_DOUBLE, i,
                 Tag::expected_solution, MPI_COMM_WORLD);
        MPI_Send(b.data(), b.size(), MPI_DOUBLE, i, Tag::b_vector,
                 MPI_COMM_WORLD);
    }

    int block_divided_len = len / blk_size;
    for (int i = 0; i < block_divided_len; i++) {
        for (int j = 0; j <= i; j++) {
            vector<double> block_ij(blk_size * blk_size, 0.0);
            // copy a block
            for (int k = 0; k < blk_size; k++) {
                for (int l = 0; l < blk_size; l++) {
                    block_ij[k * blk_size + l] =
                        full_matrix[(i * blk_size + k) * len + j * blk_size +
                                    l];
                }
            }
            int target = get_rank_by_block_id(i, j);
            MPI_Send(block_ij.data(), block_ij.size(), MPI_DOUBLE, target,
                     Tag::matrix_block, MPI_COMM_WORLD);
#ifdef DEBUG
            cout << "Sending block" << i << "," << j << " to " << target
                 << endl;
#endif
        }
    }
}

namespace std {
template <>
struct hash<pair<int, int>> {
    std::size_t operator()(const pair<int, int> &p) const noexcept {
        std::size_t h1 = std::hash<int>{}(p.first);
        std::size_t h2 = std::hash<int>{}(p.second);
        return h1 ^ (h2 << 1);  // or use boost::hash_combine
    }
};
}  // namespace std

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

    MPI_Recv(expected_solution.data(), len, MPI_DOUBLE, 0,
             Tag::expected_solution, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(b_vector.data(), len, MPI_DOUBLE, 0, Tag::b_vector, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    int block_divided_len = len / blk_size;
    for (int i = 0; i < block_divided_len; i++) {
        for (int j = 0; j <= i; j++) {
            if (get_rank_by_block_id(i, j) != rank) continue;
            vector<double> block_ij(blk_size * blk_size, 0.0);
            // copy a block
            MPI_Recv(block_ij.data(), block_ij.size(), MPI_DOUBLE, 0,
                     Tag::matrix_block, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // cout << rank << " received block " << i << "," << j << endl;
            blocks[pair<int, int>{i, j}] = move(block_ij);
        }
    }
}

void ldlt_rightlooking(const int rank, const int len, const int blk_size,
                       blocks_t &blocks) {
    int block_divided_len = len / blk_size;
    int Send_num = 0;
    double sendtime = 0;
    double recvtime = 0;
    for (int j = 0; j < block_divided_len; j++) {
        if (get_rank_by_block_id(0, j) == rank) {
            auto &Ajj = get_blocks(blocks, j, j);
            chol_block(Ajj, blk_size);  // Ajj = chol(Ajj)
            auto Ajj_inverse = trangular_inverse(Ajj, blk_size);
            for (int i = j + 1; i < block_divided_len; i++) {
                auto &Aij = get_blocks(blocks, i, j);
                auto newAij =
                    A_mul_Bt(Aij, Ajj_inverse, blk_size);  // Aij = Aij * Ajj_-t
                Aij = move(newAij);
                for (int k = j + 1; k <= i; k++) {
                    int target = get_rank_by_block_id(i, k);
                    auto &Lkj = get_blocks(blocks, k, j);
                    if (target == rank) {  // block is local
                        A_sub_B(get_blocks(blocks, i, k),
                                A_mul_Bt(Aij, Lkj, blk_size),
                                blk_size);  // Aik = Aik - Aij * Akj_t
                    } else {
                        sendtime -= MPI_Wtime();
                        MPI_Send(Aij.data(), Aij.size(), MPI_DOUBLE, target,
                                 Tag::L_matrix, MPI_COMM_WORLD);
                        MPI_Send(Lkj.data(), Lkj.size(), MPI_DOUBLE, target,
                                 Tag::L_matrix, MPI_COMM_WORLD);
                        sendtime += MPI_Wtime();
                        Send_num += 2;
                    }
                }
            }
        } else {
            for (int i = j + 1; i < block_divided_len; i++) {
                for (int k = j + 1; k <= i; k++) {
                    int target = get_rank_by_block_id(i, k);
                    if (target != rank) continue;
                    int source = get_rank_by_block_id(i, j);
                    vector<double> Aij(blk_size * blk_size, 0);
                    vector<double> Lkj(blk_size * blk_size, 0);

                    recvtime -= MPI_Wtime();
                    MPI_Recv(Aij.data(), Aij.size(), MPI_DOUBLE, source,
                             Tag::L_matrix, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(Lkj.data(), Lkj.size(), MPI_DOUBLE, source,
                             Tag::L_matrix, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    recvtime += MPI_Wtime();
                    A_sub_B(get_blocks(blocks, i, k),
                            A_mul_Bt(Aij, Lkj, blk_size),
                            blk_size);  // Aik = Aik - Aij * Akj_t
                }
            }
        }
    }
    for (int i = 0; i < world_size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (i != rank) continue;
        cout << "in rank " << i << endl;
        cout << "Number of send: " << Send_num << endl;
        cout << "send time : " << sendtime << endl;
        cout << "recv time : " << recvtime << endl;
    }
}

void verify(const vector<double> expected, const vector<double> b_vector,
            const vector<double> L) {
    const int n = expected.size();
    vector<double> solution(n, 0);
    vector<double> b(b_vector);

    for (int i = 0; i < n; i++) {
        double t = b[i];
        for (int j = 0; j < i; j++) t -= L[i * n + j] * solution[j];
        solution[i] = t / L[i * n + i];
    }

    /* backward solve D L^t x = y */
    for (int i = n - 1; i >= 0; i--) {
        double t = solution[i];

        for (int j = i + 1; j < n; j++) t -= L[j * n + i] * solution[j];

        solution[i] = t / L[i * n + i];
    }

    double e = 0;
    for (int i = 0; i < n; i++)
        e += (expected[i] - solution[i]) * (expected[i] - solution[i]);
    e = sqrt(e);
    cout << "error norm = " << e << endl;
    cout << "--- good if error is around n * 1e-16 = " << n * 1e-16
         << " or less" << endl;
}

void ldlt_leftlooking(const int rank, const int len, const int blk_size,
                      blocks_t &blocks) {
    int block_divided_len = len / blk_size;
    int Bcast_num = 0;
    double sendtime = 0;
    double recvtime = 0;
    for (int j = 0; j < block_divided_len; j++) {
        int owner = get_rank_by_block_id(j, 0);
        if (owner == rank) {  // if local is owner
            vector<MPI_Request> requests;

            auto &Ajj = get_blocks(blocks, j, j);
            for (int k = 0; k < j; k++) {
                auto &Ajk = get_blocks(blocks, j, k);
                MPI_Request r;
                sendtime -= MPI_Wtime();
                MPI_Ibcast(Ajk.data(), Ajk.size(), MPI_DOUBLE, rank,
                           MPI_COMM_WORLD, &r);
                // MPI_Bcast(Ajk.data(), Ajk.size(), MPI_DOUBLE, owner,
                // MPI_COMM_WORLD);
                sendtime += MPI_Wtime();
                requests.push_back(r);
                Bcast_num++;
                A_sub_B(Ajj, A_mul_Bt(Ajk, Ajk, blk_size),
                        blk_size);  // Ajj -= AjkAjk_t for k: 0~j-1
            }
            chol_block(Ajj, blk_size);
            auto Ajj_inv = trangular_inverse(Ajj, blk_size);
            sendtime -= MPI_Wtime();
            MPI_Request r;
            MPI_Ibcast(Ajj_inv.data(), Ajj_inv.size(), MPI_DOUBLE, rank, MPI_COMM_WORLD,
                       &r);
            // MPI_Bcast(Ajj_inv.data(), Ajj_inv.size(), MPI_DOUBLE, owner,
            // MPI_COMM_WORLD);
            sendtime += MPI_Wtime();
            requests.push_back(r);
            Bcast_num++;
            // if is local
            for (int i = j + 1; i < block_divided_len; i++) {
                if (get_rank_by_block_id(i, j) == rank) {
                    auto &Aij = get_blocks(blocks, i, j);
                    for (int k = 0; k < j; k++) {
                        auto &Aik = get_blocks(blocks, i, k);
                        auto &Ajk = get_blocks(blocks, j, k);
                        A_sub_B(Aij, A_mul_Bt(Aik, Ajk, blk_size),
                                blk_size);  // Aij -= AjkAjk_t for k: 0~j-1
                    }
                    auto newAij = A_mul_Bt(Aij, Ajj_inv, blk_size);
                    Aij = move(newAij);
                }
            }
            MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        } else {  // if not owner
            for (int k = 0; k < j; k++) {
                vector<double> Ajk(blk_size * blk_size, 0.0);
                recvtime -= MPI_Wtime();
                MPI_Request r;
                MPI_Ibcast(Ajk.data(), Ajk.size(), MPI_DOUBLE, owner,
                           MPI_COMM_WORLD, &r);
                MPI_Wait(&r, MPI_STATUSES_IGNORE);
                // MPI_Bcast(Ajk.data(), Ajk.size(), MPI_DOUBLE, owner,
                //           MPI_COMM_WORLD);
                recvtime += MPI_Wtime();
                for (int i = j + 1; i < block_divided_len; i++) {
                    if (get_rank_by_block_id(i, j) == rank) {
                        auto &Aij = get_blocks(blocks, i, j);
                        auto &Aik = get_blocks(blocks, i, k);
                        A_sub_B(Aij, A_mul_Bt(Aik, Ajk, blk_size),
                                blk_size);  // Aij -= AjkAjk_t for k: 0~j-1
                    }
                }
            }
            vector<double> Ajj(blk_size * blk_size, 0);
            recvtime -= MPI_Wtime();
            MPI_Request r;
            MPI_Ibcast(Ajj.data(), Ajj.size(), MPI_DOUBLE, owner,
                       MPI_COMM_WORLD, &r);
            MPI_Wait(&r, MPI_STATUSES_IGNORE);
            // MPI_Bcast(Ajj.data(), Ajj.size(), MPI_DOUBLE, owner,
            // MPI_COMM_WORLD);
            recvtime += MPI_Wtime();
            for (int i = j + 1; i < block_divided_len; i++) {
                if (get_rank_by_block_id(i, j) == rank) {
                    auto &Aij = get_blocks(blocks, i, j);
                    auto newAij = A_mul_Bt(Aij, Ajj, blk_size);
                    Aij = move(newAij);
                }
            }
        }
    }
    for (int i = 0; i < world_size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (i != rank) continue;
        cout << "in rank " << i << endl;
        cout << "Number of send: " << Bcast_num << endl;
        cout << "send time : " << sendtime << endl;
        cout << "recv time : " << recvtime << endl;
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

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
// parallel Cholesky decomposition
#ifdef LEFTLOOKING
    ldlt_leftlooking(rank, len, blk_size, blocks);
#else
    ldlt_rightlooking(rank, len, blk_size, blocks);
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    // mesure the time
    double duration = MPI_Wtime() - start;
    cout << duration << endl;

    // collecting result
    int block_divided_len = len / blk_size;
    if (rank == 0) {
        vector<double> full_matrix(len * len, 0.0);
        auto copy_to_full = [&](const int k, const int l,
                                const vector<double> &block) {
            for (int i = 0; i < blk_size; i++) {
                for (int j = 0; j < blk_size; j++) {
                    full_matrix[(k * blk_size + i) * len + l * blk_size + j] =
                        block[i * blk_size + j];
                }
            }
        };
        for (int i = 0; i < block_divided_len; i++) {
            for (int j = 0; j <= i; j++) {
                int source = get_rank_by_block_id(i, j);
                if (source != rank) {
                    vector<double> tmp(blk_size * blk_size, 0.0);
                    MPI_Recv(tmp.data(), tmp.size(), MPI_DOUBLE, source,
                             Tag::collect_result, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    copy_to_full(i, j, tmp);
                } else {
                    copy_to_full(i, j, get_blocks(blocks, i, j));
                }
            }
        }
        verify(expected_solution, b_vector, full_matrix);
    } else {  // non rank 0 node
        for (int i = 0; i < block_divided_len; i++) {
            for (int j = 0; j <= i; j++) {
                int source = get_rank_by_block_id(i, j);
                if (source == rank) {
                    auto &tmp = get_blocks(blocks, i, j);
                    MPI_Send(tmp.data(), tmp.size(), MPI_DOUBLE, 0,
                             Tag::collect_result, MPI_COMM_WORLD);
                }
            }
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
