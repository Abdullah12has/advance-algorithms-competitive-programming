#include <bits/stdc++.h>
using namespace std;

struct Fraction
{
    int p;
    int q;

    Fraction(int p = 0, int q = 1) : p(p), q(q)
    {
        simplify();
    }

    void simplify()
    {
        if (q == 0)
        {
            p = 0;
            q = 1;
            return;
        }

        int g = gcd(abs(p), abs(q));
        p /= g;
        q /= g;
    }

    bool operator<(const Fraction &other) const
    {
        return (long long)p * other.q < (long long)other.p * q;
    }

    bool operator==(const Fraction &other) const
    {
        return p == other.p && q == other.q;
    }

    bool operator<=(int val) const
    {
        return p <= (long long)val * q;
    }
};

class SegmentTree
{
private:
    int n;
    vector<int> tree;
    vector<int> lazy;

    void build(int node, int start, int end, const vector<int> &arr)
    {
        if (start == end)
        {
            tree[node] = arr[start];
        }
        else
        {
            int mid = (start + end) / 2;
            build(2 * node, start, mid, arr);
            build(2 * node + 1, mid + 1, end, arr);
            tree[node] = min(tree[2 * node], tree[2 * node + 1]);
        }
    }

void push(int node, int start, int end) {
    if (lazy[node] != 0) {
        // Apply lazy value to current node's minimum
        tree[node] += lazy[node]; 

        if (start != end) {
            // Pass the lazy value down to children
            lazy[2 * node] += lazy[node];
            lazy[2 * node + 1] += lazy[node];
        }
        lazy[node] = 0;
    }
}

    void updateRange(int node, int start, int end, int l, int r, int val)
    {
        push(node, start, end); // Apply pending updates before processing

        if (r < start || end < l)
        {
            return;
        }

        if (l <= start && end <= r)
        {
            tree[node] += val;
            if (start != end)
            {
                lazy[2 * node] += val;
                lazy[2 * node + 1] += val;
            }
            return;
        }

        int mid = (start + end) / 2;
        updateRange(2 * node, start, mid, l, r, val);
        updateRange(2 * node + 1, mid + 1, end, l, r, val);
        tree[node] = min(tree[2 * node], tree[2 * node + 1]);
    }

    int queryMin(int node, int start, int end, int l, int r)
    {
        if (r < start || end < l)
        {
            return INT_MAX;
        }

        push(node, start, end); // Apply pending updates before querying

        if (l <= start && end <= r)
        {
            return tree[node];
        }

        int mid = (start + end) / 2;
        int left_min = queryMin(2 * node, start, mid, l, r);
        int right_min = queryMin(2 * node + 1, mid + 1, end, l, r);
        return min(left_min, right_min);
    }

public:
    SegmentTree(const vector<int> &arr)
    {
        n = arr.size();
        tree.resize(4 * n);
        lazy.resize(4 * n, 0);
        build(1, 0, n - 1, arr);
    }

    void updateRange(int l, int r, int val)
    {
        updateRange(1, 0, n - 1, l, r, val);
    }

    int queryMin(int l, int r)
    {
        return queryMin(1, 0, n - 1, l, r);
    }

    int getMin()
    {
        push(1, 0, n - 1); // Ensure root is up-to-date
        return tree[1];
    }
};

int n;
int L;
vector<string> schedules;
int max_possible_sleep = 0;
bool max_sleep_computed = false;

vector<string> scaled_schedule;
int scaled_L;
int scaled_T;

vector<vector<pair<int, int>>> valid_intervals;
vector<pair<int, int>> sleep_assignment;

vector<pair<int, int>> FindAvailableSegments(const string &schedule)
{
    vector<pair<int, int>> segments;
    int start = -1;

    for (int i = 0; i < (int)schedule.length(); i++)
    {
        if (schedule[i] == '.')
        {
            if (start == -1)
            {
                start = i;
            }
        }
        else
        {
            if (start != -1)
            {
                segments.push_back({start, i});
                start = -1;
            }
        }
    }

    if (start != -1)
    {
        segments.push_back({start, (int)schedule.length()});
    }

    return segments;
}

int ComputeMaxPossibleSleep()
{
    if (max_sleep_computed)
    {
        return max_possible_sleep;
    }

    max_possible_sleep = INT_MAX;
    for (const string &schedule : schedules)
    {
        int best = 0;
        int current = 0;
        for (char c : schedule)
        {
            if (c == '.')
            {
                current++;
                best = max(best, current);
            }
            else
            {
                current = 0;
            }
        }
        max_possible_sleep = min(max_possible_sleep, best);
    }

    if (max_possible_sleep == INT_MAX)
    {
        max_possible_sleep = 0;
    }

    max_sleep_computed = true;
    return max_possible_sleep;
}

set<Fraction> GenerateCandidateFractions(int maxSleepTime)
{
    set<Fraction> fractions;

    if (maxSleepTime == 0)
    {
        fractions.insert(Fraction(0, 1));
        return fractions;
    }

    for (int q = 1; q <= 9; q++)
    {
        for (int p = 0; p <= maxSleepTime * q; p++)
        {
            Fraction f(p, q);
            if (f <= maxSleepTime)
            {
                fractions.insert(f);
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        auto segments = FindAvailableSegments(schedules[i]);
        for (auto [start, end] : segments)
        {
            int len = end - start;
            for (int divisor = 1; divisor <= len && divisor <= 20; divisor++)
            {
                Fraction f(len, divisor);
                if (f <= maxSleepTime)
                {
                    fractions.insert(f);
                }
            }
        }
    }

    return fractions;
}

vector<pair<int, int>> PrecomputeSleepIntervals(int caretaker_id, int T, int timeline_length)
{
    vector<pair<int, int>> intervals;

    if (T == 0)
    {
        intervals.push_back({0, 0});
        return intervals;
    }

    int start = 0;

    while (start <= timeline_length - T)
    {
        bool can_sleep = true;
        int first_blocked = -1;

        for (int t = start; t < start + T; t++)
        {
            if (scaled_schedule[caretaker_id][t] == 'X')
            {
                can_sleep = false;
                first_blocked = t;
                break;
            }
        }

        if (can_sleep)
        {
            intervals.push_back({start, start + T});
            start++;
        }
        else
        {
            start = first_blocked + 1;
        }
    }

    return intervals;
}

vector<int> InitializeCoverage(int timeline_length)
{
    vector<int> coverage(timeline_length, 0);

    for (int t = 0; t < timeline_length; t++)
    {
        for (int i = 0; i < n; i++)
        {
            if (scaled_schedule[i][t] == '.')
            {
                coverage[t]++;
            }
        }
    }

    return coverage;
}

bool Backtrack(vector<bool> &assigned, int T, int timeline_length, SegmentTree &coverage_tree)
{
    int assigned_count = 0;
    for (int i = 0; i < n; i++)
    {
        if (assigned[i])
            assigned_count++;
    }

    if (assigned_count == n)
    {
        return coverage_tree.getMin() >= 1;
    }

    if (coverage_tree.getMin() < 1)
    {
        return false;
    }

    int best_caretaker = -1;
    int min_intervals = INT_MAX;

    for (int i = 0; i < n; i++)
    {
        if (!assigned[i])
        {
            int total_intervals = valid_intervals[i].size();

            if (total_intervals == 0)
            {
                return false;
            }

            if (total_intervals < min_intervals)
            {
                min_intervals = total_intervals;
                best_caretaker = i;
            }
        }
    }

    if (best_caretaker == -1)
    {
        return false;
    }

    int caretaker_idx = best_caretaker;
    assigned[caretaker_idx] = true;

    vector<pair<int, pair<int, int>>> interval_scores;
    for (auto [start, end] : valid_intervals[caretaker_idx])
    {
        int min_coverage = INT_MAX;

        if (end > start)
        {
            min_coverage = coverage_tree.queryMin(start, end - 1);
        }

        interval_scores.push_back({-min_coverage, {start, end}});
    }

    sort(interval_scores.begin(), interval_scores.end());

    for (auto [score, interval] : interval_scores)
    {
        auto [start, end] = interval;

        sleep_assignment[caretaker_idx] = {start, end};

        if (end > start)
        {
            coverage_tree.updateRange(start, end - 1, -1);
        }

        if (coverage_tree.getMin() >= 1)
        {
            if (Backtrack(assigned, T, timeline_length, coverage_tree))
            {
                return true;
            }
        }

        if (end > start)
        {
            coverage_tree.updateRange(start, end - 1, 1);
        }
        sleep_assignment[caretaker_idx] = {-1, -1};
    }

    assigned[caretaker_idx] = false;

    return false;
}

bool IsFeasible(Fraction T)
{
    int p = T.p;
    int q = T.q;

    scaled_L = L * q;
    scaled_T = p;

    if (scaled_L == 0 && scaled_T > 0)
        return false;
    if (scaled_L == 0 && scaled_T == 0)
        return true;

    scaled_schedule.assign(n, string(scaled_L, '.'));
    for (int i = 0; i < n; i++)
    {
        for (int t = 0; t < L; t++)
        {
            for (int k = 0; k < q; k++)
            {
                scaled_schedule[i][t * q + k] = schedules[i][t];
            }
        }
    }

    valid_intervals.assign(n, vector<pair<int, int>>());
    for (int i = 0; i < n; i++)
    {
        valid_intervals[i] = PrecomputeSleepIntervals(i, scaled_T, scaled_L);

        if (valid_intervals[i].empty())
        {
            return false;
        }
    }

    vector<int> initial_coverage = InitializeCoverage(scaled_L);

    for (int cov : initial_coverage)
    {
        if (cov == 0)
            return false;
    }

    SegmentTree coverage_tree(initial_coverage);

    sleep_assignment.assign(n, {-1, -1});

    vector<bool> assigned(n, false);

    return Backtrack(assigned, scaled_T, scaled_L, coverage_tree);
}

string MaximizeSleepTime()
{
    for (int t = 0; t < L; t++)
    {
        bool all_busy = true;

        for (int i = 0; i < n; i++)
        {
            if (schedules[i][t] == '.')
            {
                all_busy = false;
                break;
            }
        }

        if (all_busy)
        {
            return "SAD CAT";
        }
    }

    int maxSleepTime = ComputeMaxPossibleSleep();

    if (maxSleepTime == 0)
    {
        Fraction zero(0, 1);
        return IsFeasible(zero) ? "0/1" : "SAD CAT";
    }

    set<Fraction> candidates = GenerateCandidateFractions(maxSleepTime);

    vector<Fraction> sorted_candidates(candidates.begin(), candidates.end());

    sort(sorted_candidates.rbegin(), sorted_candidates.rend());

    Fraction best_fraction(0, 1);
    bool found = false;
    for (const Fraction &T : sorted_candidates)
    {
        if (IsFeasible(T))
        {
            best_fraction = T;
            found = true;
            break;
        }
    }

    if (!found)
    {
        if (IsFeasible(Fraction(0, 1)))
        {
            return "0/1";
        }
        return "SAD CAT";
    }

    return to_string(best_fraction.p) + "/" + to_string(best_fraction.q);
}

int main()
{
    cin >> n >> L;
    schedules.resize(n);
    for (int i = 0; i < n; i++)
    {
        cin >> schedules[i];
    }

    cout << MaximizeSleepTime() << endl;

    return 0;
}