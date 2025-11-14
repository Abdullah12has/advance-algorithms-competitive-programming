# Dexter's Caretakers: Maximizing Sleep Time

## Problem Statement

You have **n** caretakers responsible for looking after Dexter over **L** hours (from hour 0 to hour L). Each caretaker has a schedule indicating when they are busy with chores (`X`) and when they are available (`.`). If a caretaker is not busy with chores, they can either care for Dexter or sleep.

**Goal:**
Find the **maximal sleep time T** (as a rational number in simplest form) such that:
- Every caretaker sleeps for exactly **T** hours.
- At every moment in `[0, L)`, at least one caretaker is available to care for Dexter (i.e., not busy with chores and not asleep).

If it's impossible to create such a schedule, output `SAD CAT`.

---

## Input Format

- **First line:** Two integers, `n` (number of caretakers) and `L` (number of hours).
- **Next n lines:** Each line is a string of length `L`, representing the schedule of a caretaker:
  - `X`: Busy with chores during that hour.
  - `.`: Available (can sleep or care for Dexter).

---

## Output Format

- If a valid schedule exists, output the maximal `T` as a rational number in simplest form: `p/q` (max sleep time).
- If no valid schedule exists, output `SAD CAT`.

---

## Examples

### Example 1
**Input:**
```
3 6
..X.XX
.X..X.
X..X..
```
**Output:**
```
4/3
```
**Explanation:**
Each caretaker sleeps for `4/3` hours during non-overlapping intervals, ensuring Dexter is always cared for.

---

### Example 2
**Input:**
```
3 2
..
XX
..
```
**Output:**
```
0/1
```
**Explanation:**
The second caretaker is always busy, so they cannot sleep.

---

### Example 3
**Input:**
```
1 3
.X.
```
**Output:**
```
SAD CAT
```
**Explanation:**
There is a moment when no one is available to care for Dexter.

---

## Additional Examples

| Input (n, L) | Description | Solution |
|--------------|-------------|----------|
| 6, 1000 | First caretaker is always free; others are busy except for `[123, 177)`. | `27/1` |
| 10, 30000 | No one has any chores. | `15000/1` |
| 18, 100000 | The i-th caretaker (0-indexed) is free during all intervals `[x, min(x+2, L))` where `x ≡ i mod 18`. | `2/1` |

---

## Constraints

| Subtask | Limits | Points |
|---------|--------|--------|
| 1 | `1 ≤ n ≤ 6`, `1 ≤ L ≤ 1000` | 2 |
| 2 | `1 ≤ n ≤ 10`, `1 ≤ L ≤ 30000` | 2 |
| 3 | `1 ≤ n ≤ 18`, `1 ≤ L ≤ 100000` | 6 |

---

## Key Points

- **Sleep Interval:** Each caretaker can sleep for `T` hours, starting at any time `a` such that `a + T ≤ L`.
- **Availability:** At every moment, at least one caretaker must be available (not busy with chores and not asleep).
- **Rational Output:** The answer must be a rational number in simplest form (`p/q`).
- **Impossible Case:** If no valid schedule exists, output `SAD CAT`.
