// Name:       Shuang Li
// Login ID:   shangl3
// Student ID: 1044137

#include <float.h> // for DBL_MAX
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h> // for clock()
#include <time.h>     // for clock()
#include <assert.h>

// cpp headers
#include <iostream>
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

double cell_cost(int x, int y, params* par, double** board, double** cache) {
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

    memset(cand_data, 0, y_size * x_size * sizeof(T));

    for (int i = 0; i < x_size; i++) {
        cand[i] = cand_data + i * y_size;
    }

    return cand;
}

node** init_cand(int x_size, int y_size) {
    node** cand = init_2D<node>(x_size, y_size);
    for (int i = 0; i < x_size; i++) {
        for (int j = 0; j < y_size; j++)
            cand[i][j].cost = DBL_MAX;
    }
    return cand;
}

/******************************************************************************/

/* UNSEEN must be 0, as nodes initialized using memset */
#define UNSEEN 0
#define OPEN 1
#define CLOSED 2

void a_star(double** board, int x_size, int y_size, params par) {
    int x_end = x_size - 1;
    int y_end = y_size - 1;

    priority_queue<node*, vector<node*>, Compare> pq;
    double** cell_cost_cache = init_2D<double>(x_size, y_size);

    node* pivot;
    node** cand = init_cand(x_size, y_size);
    pq.push(&(cand[0][0]));
    cand[0][0].cost = cell_cost(0, 0, &par, board, cell_cost_cache);

    while (!is_equal(pivot = pop(&pq), x_end, y_end)) {
        pivot->is_closed = CLOSED;

        /* Expand all neighbours */
#pragma omp parallel for collapse(2)
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                int new_x = pivot->x + dx;
                int new_y = pivot->y + dy;
                if (new_x < 0 || new_x > x_end || new_y < 0 || new_y > y_end ||
                    (dx == 0 && dy == 0))
                    continue;
                if (!cand[new_x][new_y].is_closed) {
                    /* Note: this calculates costs multiple times */
                    /* You will probably want to avoid that, */
                    /* but this version is easy to parallelize. */
                    double node_cost =
                        cell_cost(new_x, new_y, &par, board, cell_cost_cache);
                    if (pivot->cost + node_cost < cand[new_x][new_y].cost) {
                        cand[new_x][new_y].cost = pivot->cost + node_cost;
                        cand[new_x][new_y].x = new_x;
                        cand[new_x][new_y].y = new_y;
                        cand[new_x][new_y].dx = dx;
                        cand[new_x][new_y].dy = dy;
                        /* Here we simply insert a better path into the PQ. */
                        /* It is more efficient to change the weight of */
                        /* the old entry, but this also works. */
#pragma omp critical
                        pq.push(&(cand[new_x][new_y]));
                    }
                }
            }
        }
        DEBUG(for (int i = 0; i < y_size; i++) {
            for (int j = 0; j < x_size; j++) {
                printf("(%lg)%lg%c ", board[i][j], cand[i][j].cost,
                       " _*"[cand[i][j].is_closed]);
            }
            printf("\n");
        })
    }

    node* p = &cand[x_end][y_end];
    while (!is_equal(p, 0, 0)) {
        printf("%d %d %g %g\n", p->x, p->y, board[p->x][p->y], p->cost);
        p = &(cand[p->x - p->dx][p->y - p->dy]);
    }
    printf("%d %d %g %g\n", 0, 0, board[0][0], p->cost);
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
