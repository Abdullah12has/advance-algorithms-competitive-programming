#include <bits/stdc++.h>
using namespace std;

const int MAXN = 999999;
const long long INF = 1e18;

int n, k;
string s;

// Prefix sum array to calculate balance efficiently
// balance_prefix[i] = number of '<' minus number of '>' in s[0..i-1]
vector<int> balance_prefix;

/**
 * Precomputes the balance prefix array
 * Balance increases by 1 for '<' and decreases by 1 for '>'
 * This allows O(1) balance calculation for any substring
 */
void precompute_balance() {
    balance_prefix.resize(n + 1);
    balance_prefix[0] = 0;
    
    for (int i = 0; i < n; i++) {
        if (s[i] == '<')
            balance_prefix[i + 1] = balance_prefix[i] + 1;
        else
            balance_prefix[i + 1] = balance_prefix[i] - 1;
    }
    
    // VISUALIZATION: Print balance prefix array
    cout << "Balance Prefix Array: ";
    for (int i = 0; i <= n; i++) {
        cout << balance_prefix[i] << " ";
    }
    cout << "\n\n";
}

// Memoization map for cost function
// Key: encoded as l * MAXN + r for range [l, r]
unordered_map<int, long long> costMemo;

/**
 * Counts the number of valid balanced substrings in range [l, r]
 * A substring is valid if it's balanced (equal '<' and '>')
 * and never goes negative (no prefix has more '>' than '<')
 * 
 * This is the COST function for having a segment [l, r]
 * Cost = number of valid substrings NOT used in the partition
 */
long long count_valid(int l, int r) {
    if (l > r)
        return 0;

    // Create unique key for memoization
    int key = l * MAXN + r;
    
    // Return cached result if available
    if (costMemo.count(key))
        return costMemo[key];

    long long count = 0;
    
    // Count all balanced substrings in [l, r]
    // O(n^2) loop - tries all possible substrings
    for (int i = l; i <= r; i++) {
        int balance = 0;
        
        for (int j = i; j <= r; j++) {
            // Calculate balance for substring [i, j]
            balance = balance_prefix[j + 1] - balance_prefix[i];
            
            // Early termination: if balance goes negative, 
            // all longer substrings starting at i will also be invalid
            if (balance < 0)
                break;
            
            // If balanced (balance = 0), count this substring
            if (balance == 0)
                count++;
        }
    }

    // Cache and return result
    costMemo[key] = count;
    return count;
}

// DP arrays
// dp[i] = minimum cost to partition s[0..i] into some number of segments
long long dp[MAXN];
long long new_dp[MAXN];  // Temporary array for computing next DP layer

/**
 * Divide and Conquer optimization for DP computation
 * Uses convex hull trick / monotone optimal split point property
 * 
 * Parameters:
 * j      - current number of partitions
 * l, r   - range of positions we're computing for
 * opt_l, opt_r - range where optimal split point can be found
 * 
 * Key insight: if k* is optimal split for position i, 
 * then optimal split for i+1 is >= k*
 */
void compute(int j, int l, int r, int opt_l, int opt_r) {
    if (l > r)
        return;

    int mid = (l + r) / 2;
    long long best = INF;
    int best_k = -1;

    cout << "  Computing position " << mid << " (range [" << l << ", " << r 
         << "]), searching splits in [" << opt_l << ", " << opt_r << "]\n";

    // Try all possible split points in the optimal range
    // Split at position k means: [0..k] is handled by previous partitions,
    // [k+1..mid] is a new segment
    for (int k = opt_l; k <= min(mid - 1, opt_r); k++) {
        long long curr = dp[k] + count_valid(k + 1, mid);
        
        if (curr < best) {
            best = curr;
            best_k = k;
        }
    }

    new_dp[mid] = best;
    
    cout << "    Best split for position " << mid << " is at " << best_k 
         << " with cost " << best << "\n";

    // Recursively compute left and right halves
    // Key optimization: constrain search range using monotone property
    if (best_k != -1) {
        compute(j, l, mid - 1, opt_l, best_k);      // Left half: splits <= best_k
        compute(j, mid + 1, r, best_k, opt_r);      // Right half: splits >= best_k
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    cin >> n >> k >> s;

    cout << "=== INPUT ===" << endl;
    cout << "n = " << n << ", k = " << k << endl;
    cout << "s = " << s << "\n\n";

    // Precompute balance array for efficient range queries
    precompute_balance();

    // Initialize DP arrays with infinity
    for (int i = 0; i < n; i++) {
        dp[i] = INF;
        new_dp[i] = INF;
    }

    // BASE CASE: partition string into 1 segment
    // dp[i] = cost of having entire s[0..i] as one segment
    cout << "=== BASE CASE (j=1) ===" << endl;
    for (int i = 0; i < n; i++) {
        dp[i] = count_valid(0, i);
        if (i < 10 || i == n-1) {  // Print first 10 and last
            cout << "dp[" << i << "] = " << dp[i] << " (s[0.." << i << "])\n";
        }
    }
    cout << "\n";

    // DP LAYERS: compute for j = 2 to k partitions
    for (int j = 2; j <= k; j++) {
        cout << "=== LAYER j=" << j << " (partitioning into " << j << " segments) ===" << endl;
        
        // Reset new_dp array
        for (int i = 0; i < n; i++) {
            new_dp[i] = INF;
        }

        // Compute DP values using divide-and-conquer optimization
        // Positions [0, j-2] cannot be split into j parts, so start from j-1
        compute(j, j - 1, n - 1, j - 2, n - 2);
        
        // Copy new_dp to dp for next iteration
        for (int i = 0; i < n; i++) {
            dp[i] = new_dp[i];
        }
        
        // Print some results
        cout << "Results after layer " << j << ":\n";
        for (int i = j-1; i < min(j+5, n); i++) {
            cout << "  dp[" << i << "] = " << dp[i] << "\n";
        }
        if (n-1 >= j+5) {
            cout << "  ...\n";
            cout << "  dp[" << n-1 << "] = " << dp[n-1] << "\n";
        }
        cout << "\n";
        
        // Periodically clear memo to save memory
        if (j % 10 == 0) {
            costMemo.clear();
            cout << "  [Cleared memoization cache]\n\n";
        }
    }

    cout << "=== FINAL ANSWER ===" << endl;
    cout << "Minimum cost to partition into " << k << " segments: " << dp[n - 1] << endl;

    return 0;
}