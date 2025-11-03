#include <bits/stdc++.h>
using namespace std;

const long long INF = 1e18;
int n, k;
string s;

vector<int> balance_prefix;
vector<vector<int>> pos;
vector<int> k_j;

void precompute_balance() {
    balance_prefix.resize(n + 1);
    balance_prefix[0] = 0;
    for (int i = 0; i < n; i++) {
        balance_prefix[i + 1] = balance_prefix[i] + (s[i] == '<' ? 1 : -1);
    }
}

class MaxSegmentTree {
public:
    vector<int> tree;
    int size;

    MaxSegmentTree(int _size) {
        size = _size;
        tree.assign(4 * size, -1);
    }

    void update(int idx, int val, int node = 1, int left = 0, int right = -1) {
        if (right == -1) right = size - 1;
        if (left == right) {
            tree[node] = max(tree[node], val);
            return;
        }
        int mid = (left + right) / 2;
        if (idx <= mid) update(idx, val, node * 2, left, mid);
        else update(idx, val, node * 2 + 1, mid + 1, right);
        tree[node] = max(tree[node * 2], tree[node * 2 + 1]);
    }

    int query(int qleft, int qright, int node = 1, int left = 0, int right = -1) {
        if (right == -1) right = size - 1;
        if (qleft > right || qright < left) return -1;
        if (qleft <= left && right <= qright) return tree[node];
        int mid = (left + right) / 2;
        int p1 = query(qleft, qright, node * 2, left, mid);
        int p2 = query(qleft, qright, node * 2 + 1, mid + 1, right);
        return max(p1, p2);
    }
};

void precompute_k_j() {
    int OFFSET = n;
    pos.assign(2 * OFFSET + 1, vector<int>());
    
    for (int i = 0; i <= n; i++) {
        pos[balance_prefix[i] + OFFSET].push_back(i);
    }

    k_j.assign(n, -1);
    MaxSegmentTree st(2 * OFFSET + 1);
    st.update(balance_prefix[0] + OFFSET, 0);

    for (int m = 1; m <= n; m++) {
        int h = balance_prefix[m];
        int h_ = h + OFFSET;
        int q_r = h_ - 1;
        int q_l = 0;
        int max_pos = (q_l <= q_r) ? st.query(q_l, q_r) : -1;
        if (m - 1 >= 0) k_j[m - 1] = max_pos;
        st.update(h_, m);
    }
}

unordered_map<long long, long long> costMemo;

long long count_valid(int l, int r) {
    if (l > r) return 0;

    long long key = ((long long)l << 20) | r; 
    
    auto it = costMemo.find(key);
    if (it != costMemo.end()) return it->second;

    long long count = 0;
    int OFFSET = n;

    for (int j = l; j <= r; j++) {
        int h = balance_prefix[j + 1];
        int low = max(l, k_j[j]);
        int h_ = h + OFFSET;

        auto &v = pos[h_];
        auto it1 = lower_bound(v.begin(), v.end(), low);
        auto it2 = lower_bound(it1, v.end(), j + 1);
        count += distance(it1, it2);
    }

    costMemo[key] = count;
    return count;
}

struct DP_Result {
    long long cost;
    int count;
};

struct Line {
    long long m, b;
    mutable double x_intersect;
    
    bool operator<(const Line& o) const {
        return m < o.m;
    }
    
    bool operator<(double x) const {
        return x_intersect < x;
    }
};

struct LineContainer : multiset<Line, less<>> {
    const long long INF_COORD = LLONG_MAX;
    
    double intersect(const Line& a, const Line& b) {
        if (a.m == b.m) return a.b > b.b ? INF_COORD : -INF_COORD;
        return (double)(b.b - a.b) / (a.m - b.m);
    }
    
    bool bad(iterator y) {
        auto z = next(y);
        if (y == begin()) {
            if (z == end()) return false;
            return y->m == z->m && y->b <= z->b;
        }
        auto x = prev(y);
        if (z == end()) return y->m == x->m && y->b <= x->b;
        return intersect(*x, *z) <= intersect(*x, *y);
    }
    
    void add(long long m, long long b) {
        auto y = insert({m, b, 0});
        y->x_intersect = intersect(*y, *next(y));
        if (bad(y)) { erase(y); return; }
        while (next(y) != end() && bad(next(y))) erase(next(y));
        while (y != begin() && bad(prev(y))) erase(prev(y));
        auto x = y;
        while (x != begin()) {
            x--;
            x->x_intersect = intersect(*x, *next(x));
        }
    }
    
    long long query(long long x) {
        if (empty()) return INF;
        auto l = lower_bound((double)x);
        return l->m * x + l->b;
    }
};

DP_Result solve_dp_with_lambda_cht(long long lambda) {
    vector<long long> G(n + 1, INF);
    vector<int> C(n + 1, n + 1);
    
    G[0] = 0;
    C[0] = 0;
    
    LineContainer hull;
    hull.add(0, 0); 
    
    for (int i = 1; i <= n; i++) {
        long long best_cost = INF;
        int best_count = n + 1;
        
        
        for (int x = 0; x < i; x++) {
            if (G[x] == INF) continue;
            
            long long current_cost = G[x] + count_valid(x, i - 1) + lambda;
            
            if (current_cost < best_cost) {
                best_cost = current_cost;
                best_count = C[x] + 1;
            } else if (current_cost == best_cost) {
                best_count = min(best_count, C[x] + 1);
            }
        }
        
        G[i] = best_cost;
        C[i] = best_count;
    }
    
    return {G[n], C[n]};
}

long long solve_ctml_optimized() {
    precompute_balance();
    precompute_k_j();
    
    // Binary search on lambda
    long long lambda_low = 0;
    long long lambda_high = (long long)n * n;
    long long final_lambda = 0;
    
    DP_Result best_result = {INF, n + 1};

    while (lambda_low <= lambda_high) {
        long long lambda_mid = lambda_low + (lambda_high - lambda_low) / 2;
        
        DP_Result current_result = solve_dp_with_lambda_cht(lambda_mid);
        
        if (current_result.count <= k) {
            best_result = current_result;
            final_lambda = lambda_mid;
            lambda_high = lambda_mid - 1;
        } else {
            lambda_low = lambda_mid + 1;
        }
    }
    
    // Recalculate with final lambda to get accurate result
    DP_Result final_dp = solve_dp_with_lambda_cht(final_lambda);
    
    long long true_min_cost = final_dp.cost - (final_lambda * final_dp.count);
    
    return true_min_cost;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    
    cin >> n >> k >> s;
    
    cout << solve_ctml_optimized() << endl;
    
    return 0;
}