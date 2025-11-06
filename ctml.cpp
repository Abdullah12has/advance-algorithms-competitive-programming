#include <bits/stdc++.h>
using namespace std;

int n, k;
string s;
vector<vector<long long>> dp;
vector<vector<int>> opt;
vector<vector<long long>> cost;

void precompute_costs() {
    cost.assign(n, vector<long long>(n, 0));
    
    for (int l = 0; l < n; l++) {
        stack<long long> st;
        long long dp_val = 0;
        long long last = 0;
        
        for (int r = l; r < n; r++) {
            if (s[r] == '<') {
                st.push(last);
                last = 0;
            } else { 
                if (!st.empty()) {
                    long long l_val = st.top();
                    st.pop();
                    last = l_val + 1;
                    dp_val += last;
                } else {
                    last = 0;
                }
            }
            cost[l][r] = dp_val;
        }
    }
}

long long C(int i, int j) {
    if (i > j) return 0;
    return cost[i][j];
}

long long solve() {
    precompute_costs();
    

    dp.assign(k + 1, vector<long long>(n + 1, LLONG_MAX));
    opt.assign(k + 1, vector<int>(n + 1, 0));
    

    for (int j = 1; j <= n; j++) {
        dp[1][j] = C(0, j - 1);
        opt[1][j] = 0;
    }
    
    //  Knuth's optimization

    for (int i = 2; i <= k; i++) {
        for (int j = i; j <= n; j++) {

            int lower = (j == i) ? (i - 1) : opt[i][j - 1];
            int upper = j - 1;
            
            long long best_val = LLONG_MAX;
            int best_split = -1;

            for (int x = lower; x <= upper; x++) {
                long long val = dp[i - 1][x] + C(x, j - 1);
                if (val < best_val) {
                    best_val = val;
                    best_split = x;
                }
            }
            
            dp[i][j] = best_val;
            opt[i][j] = best_split;
        }
    }
    
    return dp[k][n];
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    
    cin >> n >> k;
    cin >> s;
    
    cout << solve() << "\n";
    
    return 0;
}