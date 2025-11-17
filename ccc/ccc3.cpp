#include <bits/stdc++.h>
using namespace std;

/*
  Correct, optimized fractional-T solver.

  Approach:
   - Build each caretaker's '.' segments as closed intervals [a,b] (b is integer but treat as real).
   - Candidate T values where feasibility can change are derived from segment endpoints:
       * segment lengths (b-a)
       * differences between endpoints (b_i - a_j, a_i - a_j, etc.)
     (We include a reasonably comprehensive set of breakpoints.)
   - For each candidate T and each midpoint between consecutive candidates:
       * For each caretaker, form their allowed-start intervals: [segL, segR - T] for each seg with length >= T.
       * Convert allowed-start intervals to t-coverage intervals: [sL, sR + T].
       * Merge per-caretaker t-coverage intervals.
       * Intersect the per-caretaker merged lists across caretakers. If the final intersection is non-empty,
         there exists some time t at which ALL caretakers *could* be asleep simultaneously -> T is INVALID.
       * So feasible(T) = not(invalid).
   - Choose the maximum feasible T among tested points.
   - Convert the resulting floating T to a reduced rational using a continued fraction approximation with a
     large maximum denominator (1e6) to get an exact-looking p/q (for contest-style inputs the exact rational
     will be found).
*/

const double EPS = 1e-12;

int n, L;
vector<string> S;
vector<vector<pair<double,double>>> segs; // per-caretaker segments [a,b) but we treat as [a,b] for reals

// Merge a vector of closed intervals (assumes input intervals may overlap and are not necessarily sorted)
// Output merged closed intervals (vector of pairs [l,r]).
vector<pair<double,double>> mergeIntervals(vector<pair<double,double>> v) {
    if (v.empty()) return {};
    sort(v.begin(), v.end());
    vector<pair<double,double>> out;
    double cl = v[0].first, cr = v[0].second;
    for (size_t i = 1; i < v.size(); ++i) {
        double l = v[i].first, r = v[i].second;
        if (l <= cr + EPS) cr = max(cr, r);
        else {
            out.emplace_back(cl, cr);
            cl = l; cr = r;
        }
    }
    out.emplace_back(cl, cr);
    return out;
}

// Intersect two lists of closed intervals (both merged & sorted). Returns merged intersection.
vector<pair<double,double>> intersectTwo(const vector<pair<double,double>>& A, const vector<pair<double,double>>& B) {
    vector<pair<double,double>> out;
    size_t i=0, j=0;
    while (i < A.size() && j < B.size()) {
        double l = max(A[i].first, B[j].first);
        double r = min(A[i].second, B[j].second);
        if (l <= r + EPS) out.emplace_back(l, r);
        if (A[i].second < B[j].second) ++i; else ++j;
    }
    return mergeIntervals(out);
}

// Check whether T is feasible: returns true if it's possible to assign starts so that
// NOT all caretakers are asleep at the same time (i.e., no time t exists where every caretaker can be asleep).
// Implementation: build per-caretaker t-coverage intervals; if intersection across caretakers is non-empty -> invalid -> return false.
bool feasible(double T) {
    if (T < 0) return false;
    // Build per-caretaker merged t-coverage intervals
    vector<vector<pair<double,double>>> tcov(n);
    for (int i = 0; i < n; ++i) {
        vector<pair<double,double>> tmp;
        for (auto &sg : segs[i]) {
            double a = sg.first;
            double b = sg.second; // seg is [a,b) originally, but treat as [a,b] continuous; using integer endpoints is fine
            double len = b - a;
            if (len + EPS < T) continue;
            double sL = a;
            double sR = b - T; // allowed start interval [sL, sR]
            // the times this caretaker can cover (be asleep) are union over starts of [s, s+T], i.e. [sL, sR+T]
            double tl = sL;
            double tr = sR + T;
            // Clip into timeline [0, L] for numerical stability
            tl = max(tl, 0.0);
            tr = min(tr, (double)L);
            if (tl <= tr + EPS) tmp.emplace_back(tl, tr);
        }
        if (tmp.empty()) {
            // This caretaker cannot sleep at length T anywhere -> must be feasible (they can't be forced asleep), but
            // careful: if any caretaker has no allowed starts then there is no way they can sleep T, which means T is infeasible
            return false;
        }
        tcov[i] = mergeIntervals(tmp);
    }

    // Intersect all caretakers' t-coverage intervals. If final intersection non-empty => exists t where all can be asleep -> INVALID
    vector<pair<double,double>> inter = tcov[0];
    for (int i = 1; i < n && !inter.empty(); ++i) {
        inter = intersectTwo(inter, tcov[i]);
    }
    // If inter is non-empty -> some time t is covered by all -> T invalid (not feasible)
    if (!inter.empty()) return false;
    return true;
}

// Continued fraction rational approximation of a double x with max denominator limit.
// Returns pair (p,q) reduced.
pair<long long,long long> approximate_rational(double x, long long max_den = 1000000LL) {
    if (x < 0) {
        auto pr = approximate_rational(-x, max_den);
        return {-pr.first, pr.second};
    }
    long long a0 = floor(x);
    if (fabs(x - a0) < 1e-15) return {a0, 1};
    double frac = x;
    vector<long long> cf; // continued fraction coeffs
    for (int iter=0; iter<80; ++iter) {
        long long a = floor(frac + 1e-16);
        cf.push_back(a);
        double r = frac - a;
        if (fabs(r) < 1e-15) break;
        frac = 1.0 / r;
    }
    // Build convergents and pick best with q <= max_den
    long long bestp=0, bestq=1;
    long long p0=1, q0=0, p1=cf[0], q1=1;
    if (q1 <= max_den) { bestp = p1; bestq = q1; }
    for (size_t i=1;i<cf.size();++i) {
        long long a = cf[i];
        long long p2 = a * p1 + p0;
        long long q2 = a * q1 + q0;
        if (q2 > max_den) break;
        bestp = p2; bestq = q2;
        p0 = p1; q0 = q1; p1 = p2; q1 = q2;
    }
    if (bestq == 0) { // fallback
        long long q = 1;
        long long p = llround(x * q);
        long long g = std::gcd(std::llabs(p), q);
        return {p/g, q/g};
    }
    long long g = std::gcd(std::llabs(bestp), bestq);
    return {bestp/g, bestq/g};
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (!(cin >> n >> L)) return 0;
    S.resize(n);
    for (int i = 0; i < n; ++i) cin >> S[i];

    // Quick impossibility check: if at any integer time everyone is busy -> "SAD CAT"
    for (int t = 0; t < L; ++t) {
        bool allX = true;
        for (int i = 0; i < n; ++i) if (S[i][t] == '.') { allX = false; break; }
        if (allX) {
            cout << "SAD CAT\n";
            return 0;
        }
    }

    // Build segments per caretaker as [start, end] where end is integer position (we use end as real)
    segs.assign(n, {});
    for (int i = 0; i < n; ++i) {
        int st = -1;
        for (int t = 0; t < L; ++t) {
            if (S[i][t] == '.') {
                if (st == -1) st = t;
            } else {
                if (st != -1) {
                    segs[i].emplace_back((double)st, (double)t); // [st, t)
                    st = -1;
                }
            }
        }
        if (st != -1) segs[i].emplace_back((double)st, (double)L);
    }

    // Build candidate T values
    vector<double> cand;
    cand.push_back(0.0);
    // segment lengths
    for (int i = 0; i < n; ++i) for (auto &sg : segs[i]) cand.push_back(sg.second - sg.first);
    // differences between endpoints (comprehensive)
    for (int i = 0; i < n; ++i) for (auto &a : segs[i])
    for (int j = 0; j < n; ++j) for (auto &b : segs[j]) {
        cand.push_back(a.second - b.first);
        cand.push_back(a.second - b.second);
        cand.push_back(a.first  - b.first);
        cand.push_back(a.first  - b.second);
    }

    sort(cand.begin(), cand.end());
    cand.erase(unique(cand.begin(), cand.end(), [](double a, double b){ return fabs(a-b) < 1e-12; }), cand.end());

    // We'll test each candidate > 0, and also midpoints between consecutive candidates
    double best = 0.0;
    for (double v : cand) {
        if (v <= 0) continue;
        if (feasible(v)) best = max(best, v);
    }
    for (size_t i = 0; i + 1 < cand.size(); ++i) {
        double a = cand[i], b = cand[i+1];
        double left = max(a, 0.0);
        double right = b;
        if (right - left <= 1e-15) continue;
        double mid = (left + right) / 2.0;
        if (feasible(mid)) best = max(best, mid);
    }

    // Also test a small epsilon above 0
    if (feasible(1e-12)) best = max(best, 0.0);

    // Convert best to reduced rational via continued fraction approx with reasonably large denominator
    auto pr = approximate_rational(best, 1000000LL);
    long long p = pr.first, q = pr.second;
    if (q == 0) { cout << "0/1\n"; return 0; }
    // Ensure sign and reduce
    long long g = std::gcd(std::llabs(p), q);
    p /= g; q /= g;
    if (q < 0) { q = -q; p = -p; }
    cout << p << "/" << q << "\n";
    return 0;
}
