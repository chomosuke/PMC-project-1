// Name:       Shuang Li
// Login ID:   shangl3
// Student ID: 1044137

#include <assert.h>
#include <cmath>
#include <float.h> // for DBL_MAX
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h> // for clock()
#include <time.h>     // for clock()

// cpp headers
#include <iostream>
#include <omp.h>
#include <queue>

using namespace std;

#define MAX_NODES 100000

#define _inline
// #define _inline inline

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

double cell_cost(int x, int y, double** board, double** cache, params* par) {
    if (cache[x][y] == 0) {
        cache[x][y] = cell_cost(board[x][y], par);
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

template <typename T> T** init_2D(int x_size, int y_size) {
    T** cand = (T**)malloc(x_size * sizeof(T*));
    T* cand_data = (T*)malloc(x_size * y_size * sizeof(T));
    assert_msg(cand != NULL && cand_data != NULL, "Could not allocate open");

    for (int i = 0; i < x_size; i++) {
        cand[i] = cand_data + i * y_size;
    }

    return cand;
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

double find_delta(int x_size, int y_size, double** board, double** cache,
                  params* par) {
    double max_cost = 0;
    x_size = min(x_size, 2);
    y_size = min(y_size, 2);

    // #pragma omp parallel for collapse(2)
    for (int i = 0; i < x_size; i++) {
        for (int j = 0; j < y_size; j++) {
            double cost = cell_cost(i, j, board, cache, par);
            // #pragma omp critical
            // use reduce
            if (max_cost < cost)
                max_cost = cost;
        }
    }

    return max_cost / 8;
}

_inline void move(int x, int y, int** bmap, int** bloc, double** distance,
                  double delta, vector<vector<pair<int, int>>>& b) {
    if (bmap[x][y] != -1) {
        b[bmap[x][y]][bloc[x][y]].first = -1;
    }
    int i = floor(distance[x][y] / delta);
    while (i >= b.size()) {
        // #pragma omp critical
        b.push_back(vector<pair<int, int>>());
    }
    // #pragma omp critical
    b[i].push_back(pair(x, y));
    bmap[x][y] = i;
}

vector<pair<int, int>> par_relax(vector<pair<int, int>>& b, int** bmap,
                                 double** distance, bool light, double delta,
                                 int x_size, int y_size, double** board,
                                 double** cache, params* par) {
    static int** visited = init_2D<int>(x_size, y_size);
    static int v_true = 1;
    v_true++;
    vector<pair<int, int>> req;
    // #pragma omp parallel for
    for (int i = 0; i < b.size(); i++) {
        int x = b[i].first, y = b[i].second;
        bmap[x][y] = -1;
        if (x != -1) {
            // #pragma omp parallel for collapse(2)
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    int new_x = x + dx, new_y = y + dy;
                    if (new_x < x_size && new_y < y_size && new_x >= 0 &&
                        new_y >= 0) {
                        double cost =
                            cell_cost(new_x, new_y, board, cache, par);
                        if ((cost <= delta) == light) {

                            double new_cost = distance[x][y] + cost;

                            // #pragma omp critical
                            {
                                if (distance[new_x][new_y] < new_cost) {
                                    distance[new_x][new_y] = new_cost;
                                    if (visited[new_x][new_y] != v_true) {
                                        visited[new_x][new_y] = v_true;
                                        req.push_back(pair(new_x, new_y));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return req;
}

void a_star(double** board, int x_size, int y_size, params par) {
    int x_end = x_size - 1;
    int y_end = y_size - 1;

    double** cache = init_2D<double>(x_size, y_size);
    memset(cache[0], 0, y_size * x_size * sizeof(double));

    double delta = find_delta(x_size, y_size, board, cache, &par);

    double** distance = init_2D<double>(x_size, y_size);
    int** bmap = init_2D<int>(x_size, y_size);
    int** bloc = init_2D<int>(x_size, y_size);
    vector<vector<pair<int, int>>> b;
    // #pragma omp parallel for collapse(2)
    for (int i = 0; i < x_size; i++) {
        for (int j = 0; j < y_size; j++) {
            distance[i][j] = DBL_MAX;
            bmap[i][j] = -1;
        }
    }

    distance[0][0] = 0;
    move(0, 0, bmap, bloc, distance, delta, b);
    int i = 0;
    while (!(distance[x_end][y_end] != DBL_MAX && bmap[x_end][y_end] == -1)) {
        vector<pair<int, int>> removed;
        while (!b[i].empty()) {
            vector<pair<int, int>> req =
                par_relax(b[i], bmap, distance, true, delta, x_size, y_size,
                          board, cache, &par);
            // need fix
            removed.reserve(removed.size() +
                            std::distance(b[i].begin(), b[i].end()));
            removed.insert(removed.end(), b[i].begin(), b[i].end());
            b[i].clear();

            // #pragma omp parallel for
            for (int i = 0; i < req.size(); i++) {
                move(req[i].first, req[i].second, bmap, bloc, distance, delta,
                     b);
            }
        }
        printf("Hi %d\n%f %d\n", i, distance[x_end][y_end], bmap[x_end][y_end]);
        vector<pair<int, int>> req =
            par_relax(removed, bmap, distance, false, delta, x_size, y_size,
                      board, cache, &par);

        // #pragma omp parallel for
        for (int i = 0; i < req.size(); i++) {
            move(req[i].first, req[i].second, bmap, bloc, distance, delta, b);
        }
        i++;
    }

    free(cache[0]);
    free(cache);
    free(distance[0]);
    free(distance);
    free(bmap[0]);
    free(bmap);
    free(bloc[0]);
    free(bloc);
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
