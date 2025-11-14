#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <numeric>
#include <string>
#include <climits>
using namespace std;

struct Fraction
{
    int p, q;

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
        if (q < 0)
        {
            p = -p;
            q = -q;
        }
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

    double toDouble() const
    {
        return (double)p / q;
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

    void push(int node)
    {
        if (lazy[node] != 0)
        {
            tree[2 * node] += lazy[node];
            tree[2 * node + 1] += lazy[node];
            lazy[2 * node] += lazy[node];
            lazy[2 * node + 1] += lazy[node];
            lazy[node] = 0;
        }
    }

    void updateRange(int node, int start, int end, int l, int r, int val)
    {
        if (r < start || end < l)
            return;
        if (l <= start && end <= r)
        {
            tree[node] += val;
            lazy[node] += val;
            return;
        }
        push(node);
        int mid = (start + end) / 2;
        updateRange(2 * node, start, mid, l, r, val);
        updateRange(2 * node + 1, mid + 1, end, l, r, val);
        tree[node] = min(tree[2 * node], tree[2 * node + 1]);
    }

    int queryMin(int node, int start, int end, int l, int r)
    {
        if (r < start || end < l)
            return INT_MAX;
        if (l <= start && end <= r)
            return tree[node];
        push(node);
        int mid = (start + end) / 2;
        return min(queryMin(2 * node, start, mid, l, r),
                   queryMin(2 * node + 1, mid + 1, end, l, r));
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
        return queryMin(0, n - 1);
    }
};

int n, L;
vector<string> schedules;
vector<string> scaled_schedule;
int scaled_L, scaled_T;
vector<vector<pair<int, int>>> valid_intervals;
vector<pair<int, int>> sleep_assignment;

vector<pair<int, int>> FindAvailableSegments(const string &schedule)
{
    vector<pair<int, int>> segments;
    int start = -1;

    for (int i = 0; i < schedule.length(); i++)
    {
        if (schedule[i] == '.')
        {
            if (start == -1)
                start = i;
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

set<Fraction> GenerateCandidateFractions(int L)
{
    set<Fraction> fractions;

    for (int q = 1; q <= L && q <= 20; q++)
    {
        for (int p = 0; p <= L * q; p++)
        {
            Fraction f(p, q);
            if (f <= L)
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
                if (f <= L)
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

    for (int start = 0; start <= timeline_length - T; start++)
    {
        bool can_sleep = true;
        for (int t = start; t < start + T; t++)
        {
            if (scaled_schedule[caretaker_id][t] == 'X')
            {
                can_sleep = false;
                break;
            }
        }

        if (can_sleep)
        {
            intervals.push_back({start, start + T});
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

bool Backtrack(int caretaker_idx, int T, int timeline_length, SegmentTree &coverage_tree)
{
    if (caretaker_idx == n)
    {
        return coverage_tree.getMin() >= 1;
    }

    if (coverage_tree.getMin() < 1)
    {
        return false;
    }

    for (auto [start, end] : valid_intervals[caretaker_idx])
    {
        sleep_assignment[caretaker_idx] = {start, end};

        if (end > start)
        {
            coverage_tree.updateRange(start, end - 1, -1);
        }

        if (coverage_tree.getMin() >= 1)
        {
            if (Backtrack(caretaker_idx + 1, T, timeline_length, coverage_tree))
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

    return false;
}

bool IsFeasible(Fraction T)
{
    int p = T.p, q = T.q;

    scaled_L = L * q;
    scaled_T = p;

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
    SegmentTree coverage_tree(initial_coverage);

    sleep_assignment.assign(n, {-1, -1});

    return Backtrack(0, scaled_T, scaled_L, coverage_tree);
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

    set<Fraction> candidates = GenerateCandidateFractions(L);
    vector<Fraction> sorted_candidates(candidates.begin(), candidates.end());
    sort(sorted_candidates.rbegin(), sorted_candidates.rend());

    bool found = false;
    Fraction best_fraction(0, 1);
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
