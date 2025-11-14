#include <bits/stdc++.h>
using namespace std;

int n, k;
string s;
vector<vector<long long>> dp;
vector<vector<int>> opt;

long long C(int a, int b) {
    if (a > b) return 0;
    vector<long long> st;
    long long last = 0, res = 0;
    for (int i = a; i <= b; i++) {
        if (s[i] == '<') st.push_back(last), last = 0;
        else if (!st.empty()) {
            long long l_val = st.back();
            st.pop_back();
            res += (last = l_val + 1);
        } else last = 0;
    }
    return res;
}

void compute_dp(int i, int l, int r, int opt_l, int opt_r) {
    if (l > r) return;
    int mid = (l + r) / 2;
    long long best_val = LLONG_MAX;
    int best_split = -1;
    for (int x = opt_l; x <= min(mid - 1, opt_r); x++) {
        long long val = dp[i - 1][x] + C(x, mid - 1);
        if (val < best_val) best_val = val, best_split = x;
    }
    dp[i][mid] = best_val;
    if (best_split != -1) {
        compute_dp(i, l, mid - 1, opt_l, best_split);
        compute_dp(i, mid + 1, r, best_split, opt_r);
    }
}

struct DPResult { long long cost; int splits; };

DPResult solve_with_lambda(long long lambda) {
    vector<long long> G(n + 1, LLONG_MAX);
    vector<int> cnt(n + 1, n + 1);
    G[0] = cnt[0] = 0;
    for (int i = 1; i <= n; i++) {
        long long best_cost = LLONG_MAX;
        int best_cnt = n + 1;
        for (int x = 0; x < i; x++) {
            if (G[x] == LLONG_MAX) continue;
            long long cost = G[x] + C(x, i - 1) + lambda;
            if (cost < best_cost || (cost == best_cost && cnt[x] + 1 < best_cnt))
                best_cost = cost, best_cnt = cnt[x] + 1;
        }
        G[i] = best_cost, cnt[i] = best_cnt;
    }
    return {G[n], cnt[n]};
}

long long solve_parametric() {
    long long low = 0, high = (long long)n * n, best_lambda = 0;
    while (low <= high) {
        long long mid = low + (high - low) / 2;
        DPResult res = solve_with_lambda(mid);
        if (res.splits <= k) best_lambda = mid, high = mid - 1;
        else low = mid + 1;
    }
    DPResult res = solve_with_lambda(best_lambda);
    return res.cost - best_lambda * k;
}

long long solve() {
    if (n > 10000 || k > 10000) return solve_parametric();
    dp.assign(k + 1, vector<long long>(n + 1, LLONG_MAX));
    for (int j = 1; j <= n; j++) dp[1][j] = C(0, j - 1);
    if (k <= 30) {
        for (int i = 2; i <= k; i++) compute_dp(i, i, n, i - 1, n - 1);
    } else {
        opt.assign(k + 1, vector<int>(n + 1, 0));
        for (int i = 2; i <= k; i++) {
            for (int j = i; j <= n; j++) {
                int lower = (j == i) ? (i - 1) : opt[i][j - 1];
                long long best_val = LLONG_MAX;
                int best_split = -1;
                for (int x = lower; x < j; x++) {
                    long long val = dp[i - 1][x] + C(x, j - 1);
                    if (val < best_val) best_val = val, best_split = x;
                }
                dp[i][j] = best_val, opt[i][j] = best_split;
            }
        }
    }
    return dp[k][n];
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    cin >> n >> k >> s;
    cout << solve() << "\n";
}