#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <numeric>
#include <string>
#include <climits>
#include <bits/stdc++.h>
using namespace std;

/*
 * ============================================================================
 * PROBLEM OVERVIEW:
 * ============================================================================
 * We have n caretakers and a timeline of L time units.
 * Each caretaker has a schedule showing when they're busy (X) or available (.)
 * We need to find the maximum sleep time T such that:
 *   1. Each caretaker can sleep for exactly T consecutive time units
 *   2. They can only sleep during their available time
 *   3. At any point in time, at least one caretaker must be awake
 *
 * The answer can be a fraction (e.g., if we scale time by a factor)
 * ============================================================================
 */

// ============================================================================
// FRACTION STRUCTURE
// ============================================================================
// Represents a fraction p/q in simplified form
// Used to represent sleep durations that might not be whole numbers
struct Fraction
{
    int p; // Numerator
    int q; // Denominator

    // Constructor: creates a fraction p/q and simplifies it
    Fraction(int p = 0, int q = 1) : p(p), q(q)
    {
        simplify();
    }

    // Simplifies the fraction to lowest terms
    // Examples: 6/8 -> 3/4,  10/5 -> 2/1
    void simplify()
    {
        // Handle division by zero
        if (q == 0)
        {
            p = 0;
            q = 1;
            return;
        }

        // Find greatest common divisor and divide both p and q by it
        int g = gcd(abs(p), abs(q));
        p /= g;
        q /= g;

        // Ensure denominator is always positive (move sign to numerator) TODO: never really need it in this case.
        // if (q < 0)
        // {
        //     p = -p;
        //     q = -q;
        // }
    }

    // Compare two fractions: is this fraction less than other?
    // Cross multiply to avoid floating point errors: p/q < a/b  iff  p*b < a*q
    bool operator<(const Fraction &other) const
    {
        return (long long)p * other.q < (long long)other.p * q;
    }

    // Check if two fractions are equal
    bool operator==(const Fraction &other) const
    {
        return p == other.p && q == other.q;
    }

    // Check if fraction is <= an integer value
    // p/q <= val  iff  p <= val*q
    bool operator<=(int val) const
    {
        return p <= (long long)val * q;
    }

    // Convert fraction to decimal (for debugging)
    // double toDouble() const
    // {
    //     return (double)p / q;
    // }
};

// ============================================================================
// SEGMENT TREE FOR RANGE MINIMUM QUERY
// ============================================================================
// Data structure that efficiently:
//   1. Updates a range of values (add a constant to all elements in range)
//   2. Queries the minimum value in a range
// Both operations run in O(log n) time
class SegmentTree
{
private:
    int n;            // Size of the array
    vector<int> tree; // Tree structure storing minimum values
    vector<int> lazy; // Lazy propagation: pending updates not yet pushed down

    // Build the segment tree from initial array
    // Node indexing: node 1 is root, left child is 2*node, right child is 2*node+1
    void build(int node, int start, int end, const vector<int> &arr)
    {
        if (start == end)
        {
            // Leaf node: store the array value
            tree[node] = arr[start];
        }
        else
        {
            // Internal node: build left and right subtrees
            int mid = (start + end) / 2;
            build(2 * node, start, mid, arr);                     // Build left subtree
            build(2 * node + 1, mid + 1, end, arr);               // Build right subtree
            tree[node] = min(tree[2 * node], tree[2 * node + 1]); // Store minimum of children
        }
    }

    // Push down lazy updates to children
    // If this node has a pending update, apply it and push to children
    void push(int node, int start, int end)
    {
        if (lazy[node] != 0)
        {
            // Apply the lazy update to current node
            tree[node] += lazy[node];

            // If not a leaf, propagate the update to children
            if (start != end)
            {
                lazy[2 * node] += lazy[node];     // Mark left child as needing update
                lazy[2 * node + 1] += lazy[node]; // Mark right child as needing update

                // Note: We update lazy values but don't immediately update tree values
                // The tree values will be updated when we visit those nodes
                int mid = (start + end) / 2;
                tree[2 * node] += lazy[node] * ((mid - start + 1));
                tree[2 * node + 1] += lazy[node] * ((end - mid));
            }

            // Clear the lazy value for current node
            lazy[node] = 0;
        }
    }

    // Update range [l, r] by adding 'val' to all elements
    // node represents the segment [start, end]
    void updateRange(int node, int start, int end, int l, int r, int val)
    {
        // No overlap: [start, end] and [l, r] don't intersect
        if (r < start || end < l)
        {
            return;
        }

        // Complete overlap: [start, end] is completely within [l, r]
        if (l <= start && end <= r)
        {
            tree[node] += val; // Update this node's value
            lazy[node] += val; // Mark as having a lazy update
            return;
        }

        // Partial overlap: push down any pending updates and recurse
        push(node, start, end);
        int mid = (start + end) / 2;
        updateRange(2 * node, start, mid, l, r, val);         // Update left child
        updateRange(2 * node + 1, mid + 1, end, l, r, val);   // Update right child
        tree[node] = min(tree[2 * node], tree[2 * node + 1]); // Recalculate minimum
    }

    // Query minimum value in range [l, r]
    // node represents the segment [start, end]
    int queryMin(int node, int start, int end, int l, int r)
    {
        // No overlap
        if (r < start || end < l)
        {
            return INT_MAX; // Return a value that won't affect the minimum
        }

        // Complete overlap: return this node's value
        if (l <= start && end <= r)
        {
            return tree[node];
        }

        // Partial overlap: push down updates and recurse
        push(node, start, end);
        int mid = (start + end) / 2;
        int left_min = queryMin(2 * node, start, mid, l, r);
        int right_min = queryMin(2 * node + 1, mid + 1, end, l, r);
        return min(left_min, right_min);
    }

public:
    // Constructor: build segment tree from initial array
    SegmentTree(const vector<int> &arr)
    {
        n = arr.size();
        tree.resize(4 * n);      // Segment tree needs at most 4n nodes
        lazy.resize(4 * n, 0);   // Initialize all lazy values to 0
        build(1, 0, n - 1, arr); // Build starting from root (node 1)
    }

    // Public interface: update range [l, r]
    void updateRange(int l, int r, int val)
    {
        updateRange(1, 0, n - 1, l, r, val);
    }

    // Public interface: query minimum in range [l, r]
    int queryMin(int l, int r)
    {
        return queryMin(1, 0, n - 1, l, r);
    }

    // Get global minimum (entire array)
    int getMin()
    {
        return tree[1]; // Root always stores minimum of entire array
    }
};

// ============================================================================
// GLOBAL VARIABLES
// ============================================================================
int n;                    // Number of caretakers
int L;                    // Timeline length
vector<string> schedules; // schedules[i] = schedule of caretaker i
                          // '.' = available, 'X' = busy

// Variables for scaled timeline (when testing fraction p/q)
vector<string> scaled_schedule; // Schedules scaled by factor q
int scaled_L;                   // Scaled timeline length = L * q
int scaled_T;                   // Scaled sleep time = p

// For backtracking
vector<vector<pair<int, int>>> valid_intervals; // valid_intervals[i] = list of valid sleep
                                                // intervals for caretaker i
vector<pair<int, int>> sleep_assignment;        // sleep_assignment[i] = assigned sleep interval
                                                // for caretaker i

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * Find all maximal continuous segments of available time (dots) in a schedule
 *
 * Example: "XX...X.X..." returns [(2,5), (6,7), (8,11)]
 *          Meaning: available from time 2-4, time 6, and time 8-10
 *
 * Returns: vector of pairs (start, end) where end is EXCLUSIVE
 */
vector<pair<int, int>> FindAvailableSegments(const string &schedule)
{
    vector<pair<int, int>> segments;
    int start = -1; // Start of current available segment (-1 if not in a segment)

    // Scan through the schedule
    for (int i = 0; i < (int)schedule.length(); i++)
    {
        if (schedule[i] == '.')
        {
            // Available time
            if (start == -1)
            {
                start = i; // Start a new segment
            }
        }
        else
        {
            // Busy time ('X')
            if (start != -1)
            {
                // End the current segment
                segments.push_back({start, i});
                start = -1;
            }
        }
    }

    // Handle case where schedule ends with available time
    if (start != -1)
    {
        segments.push_back({start, (int)schedule.length()});
    }

    return segments;
}

/**
 * Generate candidate fractions to test as potential sleep times
 *
 * We generate:
 *   1. Simple fractions p/q where q is small (to keep timeline scaling manageable)
 *   2. Fractions derived from available segment lengths
 *
 * All fractions are <= L (the timeline length)
 */
set<Fraction> GenerateCandidateFractions(int L)
{
    set<Fraction> fractions; // Use set to avoid duplicates and keep sorted

    // Generate fractions p/q where q is reasonable (denominator <= 20)
    // This limits how much we scale the timeline
    for (int q = 1; q <= L && q <= 20; q++)
    {
        for (int p = 0; p <= L * q; p++)
        {
            Fraction f(p, q);
            // Only keep fractions <= L
            if (f <= L)
            {
                fractions.insert(f);
                // cout << "Generated fraction: " << f.p << "/" << f.q << endl;
            }
        }
    }

    // Add fractions based on lengths of available segments
    // If someone has an available segment of length 'len', we might be able to
    // sleep for len/k for various divisors k
    // TODO: DOES THIS ACTUALLY HELP with Speed?
    for (int i = 0; i < n; i++)
    {
        auto segments = FindAvailableSegments(schedules[i]);
        for (auto [start, end] : segments)
        {
            // cout << "Caretaker " << i << " has available segment: (" << start << ", " << end << ")" << endl;
            int len = end - start;
            for (int divisor = 1; divisor <= len && divisor <= 20; divisor++)
            {
                Fraction f(len, divisor);
                if (f <= L)
                {
                    // cout << "Generated fraction from segment: " << f.p << "/" << f.q << endl;
                    fractions.insert(f);
                }
            }
        }
    }

    return fractions;
}

/**
 * For a given caretaker, find all possible intervals where they could sleep
 * for exactly T consecutive time units
 *
 * Parameters:
 *   caretaker_id: which caretaker we're analyzing
 *   T: sleep duration (in scaled time)
 *   timeline_length: length of scaled timeline
 *
 * Returns: list of pairs (start, end) where caretaker can sleep from start to end-1
 */
vector<pair<int, int>> PrecomputeSleepIntervals(int caretaker_id, int T, int timeline_length)
{
    vector<pair<int, int>> intervals;

    // Special case: no sleep needed
    if (T == 0)
    {
        intervals.push_back({0, 0}); // Empty interval
        return intervals;
    }

    int start = 0;

    // Try each possible starting position
    while (start <= timeline_length - T)
    {
        bool can_sleep = true;
        int first_blocked = -1;

        // Check if caretaker is available for the entire interval [start, start+T)
        for (int t = start; t < start + T; t++)
        {
            if (scaled_schedule[caretaker_id][t] == 'X')
            {
                // Blocked at time t
                can_sleep = false;
                first_blocked = t;
                break;
            }
        }

        if (can_sleep)
        {
            // This is a valid sleep interval
            intervals.push_back({start, start + T});
            start++; // Check next position
        }
        else
        {
            // Can't sleep here, jump past the blocked time
            start = first_blocked + 1;
        }
    }

    return intervals;
}

/**
 * Initialize the coverage array
 * coverage[t] = number of caretakers available at time t
 *
 * This represents how many caretakers COULD be awake at each time
 * (before we assign any sleep intervals)
 */
vector<int> InitializeCoverage(int timeline_length)
{
    vector<int> coverage(timeline_length, 0);

    // For each time slot
    for (int t = 0; t < timeline_length; t++)
    {
        // Count how many caretakers are available
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

/**
 * BACKTRACKING ALGORITHM
 *
 * Try to assign sleep intervals to all caretakers such that:
 *   1. Each caretaker sleeps for exactly T time units during their available time
 *   2. At every time point, at least one caretaker is awake
 *
 * Parameters:
 *   caretaker_idx: index of current caretaker we're assigning
 *   T: sleep duration (in scaled time)
 *   timeline_length: length of scaled timeline
 *   coverage_tree: segment tree tracking coverage (how many awake caretakers at each time)
 *
 * Returns: true if we found a valid assignment, false otherwise
 */
bool Backtrack(int caretaker_idx, int T, int timeline_length, SegmentTree &coverage_tree)
{
    // BASE CASE: all caretakers have been assigned sleep intervals
    if (caretaker_idx == n)
    {
        // Check if at every time point, at least one caretaker is awake
        // This is true if the minimum coverage across all times is >= 1
        return coverage_tree.getMin() >= 1;
    }

    // PRUNING: if current coverage is already invalid, no point continuing
    if (coverage_tree.getMin() < 1)
    {
        return false; // Some time point has 0 caretakers awake
    }

    // Try each valid sleep interval for this caretaker
    for (auto [start, end] : valid_intervals[caretaker_idx])
    {
        // ASSIGN: this caretaker sleeps during [start, end)
        sleep_assignment[caretaker_idx] = {start, end};

        // UPDATE COVERAGE: this caretaker is now asleep during [start, end-1]
        // So we decrement the coverage count for those times
        if (end > start)
        { // Skip if empty interval (T=0)
            coverage_tree.updateRange(start, end - 1, -1);
        }

        // CHECK: is coverage still valid? (minimum >= 1)
        if (coverage_tree.getMin() >= 1)
        {
            // RECURSE: try assigning sleep intervals to remaining caretakers
            if (Backtrack(caretaker_idx + 1, T, timeline_length, coverage_tree))
            {
                return true; // Found a valid complete assignment!
            }
        }

        // BACKTRACK: this assignment didn't work, undo changes
        if (end > start)
        {
            coverage_tree.updateRange(start, end - 1, 1); // Restore coverage
        }
        sleep_assignment[caretaker_idx] = {-1, -1}; // Unassign
    }

    // Tried all possible sleep intervals for this caretaker, none worked
    return false;
}

/**
 * Check if a given sleep duration T (as a fraction) is feasible
 *
 * Strategy:
 *   1. Scale the timeline by the denominator q of the fraction p/q
 *      This converts fractional sleep time into integer sleep time
 *   2. Check if we can assign sleep intervals using backtracking
 *
 * Example: T = 3/2, L = 4
 *   Scale by q=2: timeline becomes length 8, sleep time becomes 3
 *   Original: "..XX" -> Scaled: "....XXXX"
 */
bool IsFeasible(Fraction T)
{
    int p = T.p; // Numerator
    int q = T.q; // Denominator

    // Calculate scaled dimensions
    scaled_L = L * q; // New timeline length
    scaled_T = p;     // New sleep duration (now an integer!)

    // Build scaled schedules by repeating each time unit q times
    scaled_schedule.assign(n, string(scaled_L, '.'));
    for (int i = 0; i < n; i++)
    {
        for (int t = 0; t < L; t++)
        {
            // Original time t becomes times [t*q, t*q+1, ..., t*q+q-1] in scaled timeline
            for (int k = 0; k < q; k++)
            {
                scaled_schedule[i][t * q + k] = schedules[i][t];
            }
        }
    }

    // Precompute valid sleep intervals for each caretaker
    valid_intervals.assign(n, vector<pair<int, int>>());
    for (int i = 0; i < n; i++)
    {
        valid_intervals[i] = PrecomputeSleepIntervals(i, scaled_T, scaled_L);

        // If any caretaker has no valid sleep intervals, it's impossible
        if (valid_intervals[i].empty())
        {
            return false;
        }
    }

    // Initialize coverage: count available caretakers at each time
    vector<int> initial_coverage = InitializeCoverage(scaled_L);

    // Build segment tree for efficient range updates and minimum queries
    SegmentTree coverage_tree(initial_coverage);

    // Initialize sleep assignment array (will be filled by backtracking)
    sleep_assignment.assign(n, {-1, -1});

    // Try to find a valid assignment using backtracking
    return Backtrack(0, scaled_T, scaled_L, coverage_tree);
}

// ============================================================================
// MAIN ALGORITHM
// ============================================================================

/**
 * Find the maximum sleep time T such that all constraints are satisfied
 *
 * Algorithm:
 *   1. Check if problem is solvable at all (is there always at least one available caretaker?)
 *   2. Generate candidate fractions to test
 *   3. Try candidates from largest to smallest
 *   4. Return the first (largest) feasible one
 *
 * Returns: string representation of answer (e.g., "3/2" or "SAD CAT" if impossible)
 */
string MaximizeSleepTime()
{
    // ========================================================================
    // STEP 1: Initial Feasibility Check
    // ========================================================================
    // Before we do anything else, check if the problem is solvable at all
    // If there's any time where ALL caretakers are busy, it's impossible
    // because we can't guarantee someone is awake at that time

    for (int t = 0; t < L; t++)
    {
        bool all_busy = true; // Assume all are busy at time t

        // Check if at least one caretaker is available
        for (int i = 0; i < n; i++)
        {
            if (schedules[i][t] == '.')
            {
                all_busy = false; // Found an available caretaker
                break;
            }
        }

        // If all caretakers are busy at time t, problem is unsolvable
        if (all_busy)
        {
            return "SAD CAT";
        }
    }

    // ========================================================================
    // STEP 2: Generate Candidate Fractions
    // ========================================================================
    // Create a set of all fractions we want to test as potential sleep times
    set<Fraction> candidates = GenerateCandidateFractions(L);

    // Convert set to vector for easier iteration
    vector<Fraction> sorted_candidates(candidates.begin(), candidates.end());

    // Sort in DESCENDING order (largest fractions first)
    // We want to find the MAXIMUM sleep time, so we test from largest to smallest
    sort(sorted_candidates.rbegin(), sorted_candidates.rend());

    // ========================================================================
    // STEP 3: Test Candidates from Largest to Smallest
    // ========================================================================
    // Try each candidate fraction, starting with the largest
    // The first one that works is our answer (since they're sorted descending)

    bool found = false;
    Fraction best_fraction(0, 1); // Initialize to 0/1

    for (const Fraction &T : sorted_candidates)
    {
        // Test if this sleep duration T is feasible
        if (IsFeasible(T))
        {
            // Found a feasible sleep time!
            best_fraction = T;
            found = true;
            break; // Stop searching (this is the maximum)
        }
    }

    // ========================================================================
    // STEP 4: Return Result
    // ========================================================================

    if (!found)
    {
        // No feasible sleep time found
        return "SAD CAT";
    }

    // Return the fraction as a string "p/q"
    return to_string(best_fraction.p) + "/" + to_string(best_fraction.q);
}

// ============================================================================
// MAIN FUNCTION
// ============================================================================

int main()
{
    // Read input
    cin >> n >> L; // n = number of caretakers, L = timeline length
    schedules.resize(n);
    // Read each caretaker's schedule
    for (int i = 0; i < n; i++)
    {
        cin >> schedules[i];
        // schedules[i] is a string of length L
        // '.' means available, 'X' means busy
    }

    //print all the schedules
    // cout << "Schedules:" << endl;
    // for (int i = 0; i < n; i++)
    // {
    //     cout << "Caretaker " << i << ": " << schedules[i] << endl;
    // }

    // Solve the problem and output the result
    cout << MaximizeSleepTime() << endl;

    return 0;


}
