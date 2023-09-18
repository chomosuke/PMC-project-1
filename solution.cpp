// Name:       Shuang Li
// Login ID:   shangl3
// Student ID: 1044137

#include <assert.h>
#include <cmath>
#include <float.h> // for DBL_MAX
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h> // for clock()
#include <time.h>     // for clock()

// cpp headers
#include <iostream>
#include <mutex>
#include <omp.h>
#include <queue>
#include <shared_mutex>
#include <vector>

using namespace std;

#define MAX_NODES 100000

// #define _inline
#define _inline inline

// #define DEBUG(x) x
#define DEBUG(x)

// From
// https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows
double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void assert_msg(int cond, const char* msg) {
    if (!cond) {
        fprintf(stderr, "%s\n", msg);
        exit(1);
    }
}

typedef struct {
    int x, y;
    char dx, dy;
    char is_closed;
    double cost;
} node;

_inline int is_equal(node* n, int x, int y) {
    // Return true if n == NULL to ensure while loop terminates
    // even if pq_pop_min is modified to return NULL for an empty queue.
    return !n || (n->x == x && n->y == y);
}

template <typename T> T** init_2D(int x_size, int y_size) {
    T** cand = (T**)malloc(x_size * sizeof(T*));
    T* cand_data = (T*)malloc(x_size * y_size * sizeof(T));
    assert_msg(cand != NULL && cand_data != NULL, "Could not allocate open");

    for (int i = 0; i < x_size; i++) {
        cand[i] = cand_data + i * y_size;
    }

    return cand;
}

/******************************************************************************/
/*
 * Calculate the cost of each square in the grid, given its seed.
 * This is deliberately expensive so that overall program run-time is not
 * dominated by overheads.
 * More computationally expensive if res is smaller.
 * Wider range of costs if scale is larger.
 *
 * Based on Park and Miller's Oct 1988 CACM random number generator
 */

typedef struct {
    int par1, par2;
} params;

double cell_cost(long int seed, params* par) {
    const unsigned long a = 16807;
    const unsigned long m = 2147483647;

    /* For debugging only */
    // return (seed);

    /* Real code */
    seed = -seed; // Make high bits non-zero
    int res = par->par1;
    int scale = par->par2;

    int cost;

    for (cost = 0; seed >> res != 0; cost++) {
        seed = (a * seed) % m;
    }

    return (10 + (cost >> (8 * sizeof(unsigned long) - res - scale))) / 10.0;
}

double** board;
int x_size;
int y_size;
params par;

double** init_cache() {
    double** cache = init_2D<double>(x_size, y_size);
    memset(cache[0], 0, y_size * x_size * sizeof(double));
    return cache;
}

double cell_cost(int x, int y) {
    static double** cache = init_cache();
    if (cache[x][y] == 0) {
        cache[x][y] = cell_cost(board[x][y], &par);
    }
    return cache[x][y];
}

/******************************************************************************/
/* Priority queue */
/* Entries are of type *node. */
/* Ordering is specified by function  greater(node *n1, node *n2) */
/* that returns 1 if *n1 < *n2 and 0 otherwise. */

class Compare {
  public:
    bool operator()(node* n1, node* n2) { return (n1->cost > n2->cost); }
};

typedef priority_queue<node*, vector<node*>, Compare> node_priority_queue;

node* pop(node_priority_queue* pq) {
    node* t = pq->top();
    pq->pop();
    return t;
}

/******************************************************************************/

double** read_board(int x_size, int y_size) {
    double** board = (double**)malloc(x_size * sizeof(*board));
    double* board_data = (double*)malloc(x_size * y_size * sizeof(*board_data));
    assert_msg(board != NULL && board_data != NULL, "Could not allocate board");

    for (int i = 0; i < x_size; i++) {
        board[i] = board_data + i * y_size;

        for (int j = 0; j < y_size; j++) {
            assert_msg(scanf("%lf", &(board[i][j])) == 1,
                       "Failed to read board");
        }
    }

    return board;
}

node** init_cand(int x_size, int y_size) {
    node** cand = init_2D<node>(x_size, y_size);
    memset(cand[0], 0, y_size * x_size * sizeof(node));

    for (int i = 0; i < x_size; i++) {
        for (int j = 0; j < y_size; j++)
            cand[i][j].cost = DBL_MAX;
    }
    return cand;
}

/******************************************************************************/

typedef struct {
    int** in_set;
    int in_true;
    vector<pair<int, int>> set;
} set;

_inline set set_new() {
    set self;
    self.in_set = init_2D<int>(x_size, y_size);
    memset(self.in_set[0], 0, y_size * x_size * sizeof(int));
    self.in_true = 1;
    return self;
}

_inline void set_clear(set* self) {
    self->in_true++;
    self->set.clear();
}

_inline void set_insert(set* self, pair<int, int> e) {
    int x = e.first, y = e.second;
    // TODO: further examine this to see if it can be sped up
    if (self->in_set[x][y] != self->in_true) {
        self->in_set[x][y] = self->in_true;
        self->set.push_back(e);
    }
}

_inline void free_set(set self) { free(self.in_set); }

/******************************************************************************/

double find_delta() {
    double max_cost = 0;
    // int x_s = min(x_size, 2);
    // int y_s = min(y_size, 2);
    int x_s = x_size;
    int y_s = y_size;

#pragma omp parallel for collapse(2) reduction(max : max_cost)
    for (int i = 0; i < x_s; i++) {
        for (int j = 0; j < y_s; j++) {
            double cost = cell_cost(i, j);
            // use reduce
            if (max_cost < cost)
                max_cost = cost;
        }
    }

    return max_cost / 8;
}

// bmap[x][y] is the bucket (x, y) is in. -1 means not in any bucket.
int** bmap;
// This is used to fast remove (x, y) from bucket by setting b[i][bloc[x][y]] to
// -1.
int** bloc;
// when first == -1, it's a null member.
vector<vector<pair<int, int>>> b;

// must obtain b_all_lock before b_lock
shared_mutex b_all_lock;
vector<omp_lock_t> b_lock;

// if not in any bucket then it's no-op
_inline void remove_b(int x, int y) {
    if (bmap[x][y] != -1) {
        // make this member null
        b[bmap[x][y]][bloc[x][y]].first = -1;
        bmap[x][y] = -1;
    }
}

double delta;
double** dss;

_inline void insert_bi(int x, int y, int bi) {
    assert(bmap[x][y] == -1);
    bmap[x][y] = bi;
    if (bi >= b.size()) {
        b_all_lock.lock();
        do {
            vector<pair<int, int>> v;
            b.push_back(v);
            omp_lock_t l;
            omp_init_lock(&l);
            b_lock.push_back(l);
        } while (bi >= b.size());
        b_all_lock.unlock();
    }

    b_all_lock.lock_shared();
    omp_set_lock(&b_lock[bi]);
    bloc[x][y] = b[bi].size();
    b[bi].push_back(pair(x, y));
    omp_unset_lock(&b_lock[bi]);
    b_all_lock.unlock_shared();
}

// removed == NULL means heavy
void relax(vector<pair<int, int>>& bi, set* relaxed, set* removed) {
#pragma omp parallel for
    for (int j = 0; j < bi.size(); j++) {
        int x = bi[j].first, y = bi[j].second;
        if (x == -1) {
            continue;
        }

        if (removed != NULL) {
            remove_b(x, y);
#pragma omp critical(removed)
            set_insert(removed, pair(x, y));
        }

        // find neighbor where they became closer to source.
        // #pragma omp parallel for collapse(2)
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                int nx = x + dx, ny = y + dy;
                if (nx >= x_size || ny >= y_size || ny < 0 || nx < 0 ||
                    (dx == 0 && dy == 0)) {
                    continue;
                }

                // atomic read, so data read isn't corrupted, doesn't actually
                // matter if dxy has been updated by another thread or not.
                // TODO: are those atomic actually neccessary?
                // The answer is probably not, but they make the code slightly
                // more portable.
                double p_d;
#pragma omp atomic read
                p_d = dss[x][y];
                double cc = cell_cost(nx, ny);
#pragma omp critical(relaxed)
                {
                    if ((cc <= delta) == (removed != NULL)) {
                        double n_d = p_d + cc;
                        double o_d;
#pragma omp atomic read
                        o_d = dss[nx][ny];
                        if (n_d < o_d) {
#pragma omp atomic read
                            dss[nx][ny] = n_d;
                            set_insert(relaxed, pair(nx, ny));
                        }
                    }
                }
            }
        }
    }
}

void a_star(double** board_, int x_size_, int y_size_, params par_) {
    board = board_;
    x_size = x_size_;
    y_size = y_size_;
    par = par_;

    int x_end = x_size - 1;
    int y_end = y_size - 1;

    delta = find_delta();

    dss = init_2D<double>(x_size, y_size);
    bmap = init_2D<int>(x_size, y_size);
    bloc = init_2D<int>(x_size, y_size);
#pragma omp parallel for collapse(2)
    for (int i = 0; i < x_size; i++) {
        for (int j = 0; j < y_size; j++) {
            dss[i][j] = DBL_MAX;
            bmap[i][j] = -1;
        }
    }

    dss[0][0] = cell_cost(0, 0);
    insert_bi(0, 0, floor(dss[0][0] / delta));

    // loop needs to terminate once the destination has its distance figured
    // out. dest's distance will be fixed once it enters a bucket and then leave
    // that bucket. set dest_i to the bucket dest enter.
    int dest_i = INT_MAX;

    // this is set with fast clear operation
    set removed = set_new();
    set relaxed = set_new();

    for (int i = 0; i <= dest_i; i++) {
        assert(i < b.size());
        DEBUG(for (int l = 0; l < b.size(); l++) {
            cout << l << ": ";
            for (int k = 0; k < b[l].size(); k++) {
                cout << b[l][k].first << b[l][k].second;
            }
            cout << endl;
        } printf("i: %d\n\n", i);)
        while (!b[i].empty()) {
            // process b[i]
            relax(b[i], &relaxed, &removed);

            // no need to lock as this is not parallelized
            // all of b[i] should be null member now
            b[i].clear();

// now insert the relaxed ones back
#pragma omp parallel for
            for (int k = 0; k < relaxed.set.size(); k++) {
                int x = relaxed.set[k].first, y = relaxed.set[k].second;
                int bi = floor(dss[x][y] / delta);
                if (x_end == x && y_end == y) {
                    dest_i = bi;
                }
                insert_bi(x, y, bi);
            }

            set_clear(&relaxed);
        }
        // now relax the removed ones heavily
        relax(removed.set, &relaxed, NULL);
#pragma omp parallel for
        for (int k = 0; k < relaxed.set.size(); k++) {
            int x = relaxed.set[k].first, y = relaxed.set[k].second;
            int bi = floor(dss[x][y] / delta);
            if (x_end == x && y_end == y) {
                dest_i = bi;
            }
            insert_bi(x, y, bi);
        }
        set_clear(&relaxed);
        set_clear(&removed);
    }

    DEBUG(for (int x = 0; x < x_size; x++) {
        for (int y = 0; y < y_size; y++) {
            cout << dss[x][y] << ' ';
        }
        cout << endl;
    })

    // now construct the path
    int x = x_end, y = y_end;
    while (x != 0 || y != 0) {
        printf("%d %d %g %g\n", x + 1, y + 1, board[x][y], dss[x][y]);
        int n_x = x, n_y = y;
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                int nx = x + dx, ny = y + dy;
                if (nx >= x_size || ny >= y_size || ny < 0 || nx < 0 ||
                    (dx == 0 && dy == 0)) {
                    continue;
                }
                if (dss[nx][ny] < dss[n_x][n_y]) {
                    n_x = nx;
                    n_y = ny;
                }
            }
        }
        x = n_x;
        y = n_y;
    }
    printf("%d %d %g %g\n", x + 1, y + 1, board[x][y], dss[x][y]);
}

/******************************************************************************/

int main() {
    int x_size, y_size;
    double** board;
    node** open;
    int i, j;
    params par;

    assert_msg(
        scanf("%d %d %d %d", &x_size, &y_size, &(par.par1), &(par.par2)) == 4,
        "Failed to read size");
    board = read_board(x_size, y_size);

    clock_t t = clock();
    double w = get_wall_time();
    a_star(board, x_size, y_size, par);
    printf("Time: %ld %lf\n", clock() - t, get_wall_time() - w);

    return 0;
}
