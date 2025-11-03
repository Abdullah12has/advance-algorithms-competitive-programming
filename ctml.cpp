#include <bits/stdc++.h>
using namespace std;

const int MAXN = 999999;
const long long INF = 1e18;

int n, k;
string s;

// Precompute arrays 
vector<int> balance_prefix;

void precompute_balance() {
    balance_prefix.resize(n + 1);
    balance_prefix[0] = 0;
    
    for (int i = 0; i < n; i++) {
        if (s[i] == '<')
            balance_prefix[i + 1] = balance_prefix[i] + 1;
        else
            balance_prefix[i + 1] = balance_prefix[i] - 1;
    }
}

unordered_map<int, long long> costMemo; //  O(1) lookup

// need to further optimize this
long long count_valid(int l, int r) {
    if (l > r)
        return 0;

    int key = l * MAXN + r;
    
    if (costMemo.count(key))
        return costMemo[key];

    long long count = 0;
    
    // n2 -> try to optimize this somehow
    for (int i = l; i <= r; i++) {
        int balance = 0;
        
        for (int j = i; j <= r; j++) {
            balance = balance_prefix[j + 1] - balance_prefix[i];
            
            if (balance < 0)
                break; // optimization early terminztion
            
            if (balance == 0)
                count++;
        }
    }

    costMemo[key] = count;
    return count;
}

long long dp[MAXN];
long long new_dp[MAXN];

// Divide and Conquer
void compute(int j, int l, int r, int opt_l, int opt_r) {
    if (l > r)
        return;

    int mid = (l + r) / 2;
    long long best = INF;
    int best_k = -1;

    // all split points in the optimal range kinda n3 here calling count_valid
    for (int k = opt_l; k <= min(mid - 1, opt_r); k++) {
        long long curr = dp[k] + count_valid(k + 1, mid);
        if (curr < best) {
            best = curr;
            best_k = k;
        }
    }

    new_dp[mid] = best;

    if (best_k != -1) {
        compute(j, l, mid - 1, opt_l, best_k);
        compute(j, mid + 1, r, best_k, opt_r);
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    cin >> n >> k >> s;

    precompute_balance();

    // init DP 
    for (int i = 0; i < n; i++) {
        dp[i] = INF;
        new_dp[i] = INF;
    }

    // Base case
    for (int i = 0; i < n; i++) {
        dp[i] = count_valid(0, i);
    }

    // DP for j = 2 to k
    for (int j = 2; j <= k; j++) {
        for (int i = 0; i < n; i++) {
            new_dp[i] = INF;
        }

        // split at positions j-1 or later
        compute(j, j - 1, n - 1, j - 2, n - 2);
        for (int i = 0; i < n; i++) {
            dp[i] = new_dp[i];
        }
        
        if (j % 10 == 0) costMemo.clear();
    }

    cout << dp[n - 1] << endl;

    return 0;
}
