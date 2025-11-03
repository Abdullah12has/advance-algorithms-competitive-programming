#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, k;
    cin >> n >> k;

    string s;
    cin >> s;

    vector<int> match(n, -1);
    stack<int> st;
    vector<pair<int, int>> pairs;

    // Step 1: Find all valid pairs using stack
    for (int i = 0; i < n; i++) {
        if (s[i] == '<') {
            st.push(i);
        } else if (s[i] == '>' && !st.empty()) {
            int open = st.top(); st.pop();
            match[open] = i;
            pairs.push_back({open, i});
        }
    }

    // Step 2: Count how many pairs cross each position
    vector<int> crossings(n + 1, 0);
    for (auto [l, r] : pairs) {
        // Every position between l and r is crossed by this pair
        crossings[l + 1]++;   // mark start
        crossings[r]--;       // mark end
    }

    // Prefix sum to get crossing counts at each position
    for (int i = 1; i < n; i++) {
        crossings[i] += crossings[i - 1];
    }

    // Step 3: Each pair is initially valid
    int total_valid = pairs.size();

    // Step 4: Pick (k - 1) most destructive cut positions
    vector<int> cut_strength;
    for (int i = 1; i < n; i++) {
        if (crossings[i] > 0) cut_strength.push_back(crossings[i]);
    }

    sort(cut_strength.rbegin(), cut_strength.rend());

    int destroyed = 0;
    for (int i = 0; i < min(k - 1, (int)cut_strength.size()); i++) {
        destroyed += cut_strength[i];
    }

    int minimal_remaining = max(0, total_valid - destroyed);
    cout << minimal_remaining << "\n";

    return 0;
}