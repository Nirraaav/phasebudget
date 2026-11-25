// Native coverage-parameterized dynamic program for wMEC and its phase margins.
//
// This is `phasebudget.dp.CoverageDP` in C++, and it answers to it: the same
// states, the same arithmetic in the same order, and the same tie-breaking, so
// a solve here returns the same cost, haplotype, assignment and margins as the
// Python pass, bit for bit. Build with -ffp-contract=off so that no multiply
// and add are fused behind the compiler's back.
//
// Memory is O(sqrt(n) 2^k) rather than O(n 2^k): the forward pass keeps a
// table every B columns, and one descending sweep recomputes each segment of
// B tables from its checkpoint to run the traceback, the backward pass and the
// boundary costs together.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace {

using std::vector;
const double INF = std::numeric_limits<double>::infinity();

struct Problem {
    int n;
    const int32_t* offsets;  // n + 1 entries into `active`, `c0`, `c1`
    const int32_t* active;   // fragment index of each slot, column by column
    const double* c0;        // cost of giving the slot's read the value 0
    const double* c1;        // cost of giving it the value 1
    vector<int32_t> slot_in_right;

    int width(int j) const { return offsets[j + 1] - offsets[j]; }
};

// Python: `phasebudget.dp.column_sum`, Neumaier's compensated sum as CPython
// 3.12 computes `sum` over floats.
double column_sum(const double* values, int count) {
    double total = 0.0, compensation = 0.0;
    for (int i = 0; i < count; ++i) {
        const double value = values[i];
        const double partial = total + value;
        if (std::fabs(total) >= std::fabs(value))
            compensation += (total - partial) + value;
        else
            compensation += (value - partial) + total;
        total = partial;
    }
    if (compensation != 0.0 && std::isfinite(compensation)) total += compensation;
    return total;
}

// Python: `CoverageDP._column_tables` followed by `_state_costs`.
void state_costs(const Problem& p, int j, vector<double>& table) {
    const int k = p.width(j);
    const int base = p.offsets[j];
    const size_t size = size_t(1) << k;
    vector<double> delta(size, 0.0);
    const double s0 = column_sum(p.c0 + base, k);
    const double s1 = column_sum(p.c1 + base, k);
    for (int i = 0; i < k; ++i) {
        const double step = p.c1[base + i] - p.c0[base + i];
        const size_t bit = size_t(1) << i;
        for (size_t mask = bit; mask < size; ++mask)
            if (mask & bit) delta[mask] = delta[mask ^ bit] + step;
    }
    table.assign(2 * size, 0.0);
    for (size_t mask = 0; mask < size; ++mask) {
        table[mask * 2] = s0 + delta[mask];
        table[mask * 2 + 1] = s1 - delta[mask];
    }
}

// Python: `CoverageDP._shared`, as (left bits, right bits) in left-slot order.
void shared(Problem& p, int left, int right, vector<int>& lbits, vector<int>& rbits) {
    lbits.clear();
    rbits.clear();
    for (int i = 0; i < p.width(right); ++i) p.slot_in_right[p.active[p.offsets[right] + i]] = i;
    for (int i = 0; i < p.width(left); ++i) {
        const int slot = p.slot_in_right[p.active[p.offsets[left] + i]];
        if (slot >= 0) {
            lbits.push_back(i);
            rbits.push_back(slot);
        }
    }
    for (int i = 0; i < p.width(right); ++i) p.slot_in_right[p.active[p.offsets[right] + i]] = -1;
}

// Python: `CoverageDP._project` tabulated for every mask of `k` bits.
void projection(int k, const vector<int>& bits, vector<uint32_t>& code) {
    const size_t size = size_t(1) << k;
    vector<uint32_t> slot_bit(k, 0);
    for (size_t t = 0; t < bits.size(); ++t) slot_bit[bits[t]] = uint32_t(1) << t;
    code.assign(size, 0);
    for (size_t mask = 1; mask < size; ++mask)
        code[mask] = code[mask & (mask - 1)] | slot_bit[__builtin_ctzll(mask)];
}

struct Scratch {
    vector<int> lbits, rbits;
    vector<uint32_t> pl, pr;
    vector<double> best, costs;
    vector<int64_t> winner;
};

// Python: one iteration of the loop in `_run_forward`, without parents.
void forward_step(Problem& p, Scratch& s, const vector<double>& prev, int j, vector<double>& out) {
    shared(p, j - 1, j, s.lbits, s.rbits);
    projection(p.width(j - 1), s.lbits, s.pl);
    projection(p.width(j), s.rbits, s.pr);
    const size_t codes = size_t(1) << s.lbits.size();
    s.best.assign(2 * codes, INF);
    for (size_t state = 0; state < prev.size(); ++state) {
        const double cost = prev[state];
        if (cost == INF) continue;
        const size_t slot = size_t(s.pl[state >> 1]) * 2 + (state & 1);
        if (cost < s.best[slot]) s.best[slot] = cost;
    }
    state_costs(p, j, s.costs);
    out.assign(s.costs.size(), INF);
    for (size_t state = 0; state < s.costs.size(); ++state) {
        const size_t code = s.pr[state >> 1];
        const double b0 = s.best[code * 2], b1 = s.best[code * 2 + 1];
        const double b = b0 <= b1 ? b0 : b1;
        if (b == INF) continue;
        out[state] = s.costs[state] + b;
    }
}

// Python: `back[state]` of `_run_forward`, recomputed for one state only.
int64_t parent_of(Problem& p, Scratch& s, const vector<double>& prev, int j, int64_t state) {
    shared(p, j - 1, j, s.lbits, s.rbits);
    projection(p.width(j - 1), s.lbits, s.pl);
    projection(p.width(j), s.rbits, s.pr);
    const uint32_t code = s.pr[size_t(state) >> 1];
    double best[2] = {INF, INF};
    int64_t winner[2] = {-1, -1};
    for (size_t prior = 0; prior < prev.size(); ++prior) {
        const double cost = prev[prior];
        if (cost == INF || s.pl[prior >> 1] != code) continue;
        const int allele = int(prior & 1);
        if (cost < best[allele]) {
            best[allele] = cost;
            winner[allele] = int64_t(prior);
        }
    }
    return winner[best[0] <= best[1] ? 0 : 1];
}

// Python: one iteration of the loop in `_run_backward`.
void backward_step(Problem& p, Scratch& s, const vector<double>& next, int j, vector<double>& out) {
    shared(p, j, j + 1, s.lbits, s.rbits);
    projection(p.width(j), s.lbits, s.pl);
    projection(p.width(j + 1), s.rbits, s.pr);
    const size_t codes = size_t(1) << s.lbits.size();
    s.best.assign(codes, INF);
    for (size_t state = 0; state < next.size(); ++state) {
        const double cost = next[state];
        if (cost == INF) continue;
        const size_t code = s.pr[state >> 1];
        if (cost < s.best[code]) s.best[code] = cost;
    }
    state_costs(p, j, s.costs);
    out.assign(s.costs.size(), INF);
    for (size_t state = 0; state < s.costs.size(); ++state) {
        const double tail = s.best[s.pl[state >> 1]];
        if (tail < INF) out[state] = s.costs[state] + tail;
    }
}

// Python: `CoverageDP.boundary_costs`.
void boundary(Problem& p, Scratch& s, const vector<double>& fwd, const vector<double>& bwd,
              int b, double& zero, double& one) {
    shared(p, b, b + 1, s.lbits, s.rbits);
    projection(p.width(b), s.lbits, s.pl);
    projection(p.width(b + 1), s.rbits, s.pr);
    const size_t codes = size_t(1) << s.lbits.size();
    vector<double> left(2 * codes, INF), right(2 * codes, INF);
    for (size_t state = 0; state < fwd.size(); ++state) {
        const double cost = fwd[state];
        if (cost == INF) continue;
        const size_t slot = size_t(s.pl[state >> 1]) * 2 + (state & 1);
        if (cost < left[slot]) left[slot] = cost;
    }
    for (size_t state = 0; state < bwd.size(); ++state) {
        const double cost = bwd[state];
        if (cost == INF) continue;
        const size_t slot = size_t(s.pr[state >> 1]) * 2 + (state & 1);
        if (cost < right[slot]) right[slot] = cost;
    }
    double best[2] = {INF, INF};
    for (size_t code = 0; code < codes; ++code)
        for (int allele = 0; allele < 2; ++allele) {
            const double head = left[code * 2 + allele];
            if (head == INF) continue;
            for (int relation = 0; relation < 2; ++relation) {
                const double other = right[code * 2 + (allele ^ relation)];
                if (other == INF) continue;
                const double total = head + other;
                if (total < best[relation]) best[relation] = total;
            }
        }
    zero = best[0];
    one = best[1];
}

// ---- sum-product: posterior margins ---------------------------------------
//
// Python: `phasebudget.posterior.boundary_log_sums`. The same states and the
// same sweep as above, with (logsumexp, +) over minus the cost in place of
// (min, +). Sums are accumulated online as (max, scaled sum) pairs, so the
// result agrees with the Python pass to rounding, not bit for bit.

const double NEG_INF = -std::numeric_limits<double>::infinity();

struct LogSum {
    double top = NEG_INF;
    double scaled = 0.0;
    void add(double value) {
        if (value == NEG_INF) return;
        if (value > top) {
            scaled = (top == NEG_INF ? 0.0 : scaled * std::exp(top - value)) + 1.0;
            top = value;
        } else {
            scaled += std::exp(value - top);
        }
    }
    double value() const { return top == NEG_INF ? NEG_INF : top + std::log(scaled); }
};

void sum_forward_step(Problem& p, Scratch& s, const vector<double>& prev, int j, vector<double>& out) {
    shared(p, j - 1, j, s.lbits, s.rbits);
    projection(p.width(j - 1), s.lbits, s.pl);
    projection(p.width(j), s.rbits, s.pr);
    vector<LogSum> groups(size_t(1) << s.lbits.size());
    for (size_t state = 0; state < prev.size(); ++state) groups[s.pl[state >> 1]].add(prev[state]);
    state_costs(p, j, s.costs);
    out.assign(s.costs.size(), NEG_INF);
    for (size_t state = 0; state < s.costs.size(); ++state) {
        const double incoming = groups[s.pr[state >> 1]].value();
        if (incoming != NEG_INF) out[state] = -s.costs[state] + incoming;
    }
}

void sum_backward_step(Problem& p, Scratch& s, const vector<double>& next, int j, vector<double>& out) {
    shared(p, j, j + 1, s.lbits, s.rbits);
    projection(p.width(j), s.lbits, s.pl);
    projection(p.width(j + 1), s.rbits, s.pr);
    vector<LogSum> groups(size_t(1) << s.lbits.size());
    for (size_t state = 0; state < next.size(); ++state) groups[s.pr[state >> 1]].add(next[state]);
    state_costs(p, j, s.costs);
    out.assign(s.costs.size(), NEG_INF);
    for (size_t state = 0; state < s.costs.size(); ++state) {
        const double outgoing = groups[s.pl[state >> 1]].value();
        if (outgoing != NEG_INF) out[state] = -s.costs[state] + outgoing;
    }
}

void sum_boundary(Problem& p, Scratch& s, const vector<double>& fwd, const vector<double>& bwd,
                  int b, double& zero, double& one) {
    shared(p, b, b + 1, s.lbits, s.rbits);
    projection(p.width(b), s.lbits, s.pl);
    projection(p.width(b + 1), s.rbits, s.pr);
    const size_t codes = size_t(1) << s.lbits.size();
    vector<LogSum> left(2 * codes), right(2 * codes);
    for (size_t state = 0; state < fwd.size(); ++state)
        left[size_t(s.pl[state >> 1]) * 2 + (state & 1)].add(fwd[state]);
    for (size_t state = 0; state < bwd.size(); ++state)
        right[size_t(s.pr[state >> 1]) * 2 + (state & 1)].add(bwd[state]);
    LogSum total[2];
    for (size_t code = 0; code < codes; ++code)
        for (int allele = 0; allele < 2; ++allele) {
            const double head = left[code * 2 + allele].value();
            if (head == NEG_INF) continue;
            for (int relation = 0; relation < 2; ++relation) {
                const double tail = right[code * 2 + (allele ^ relation)].value();
                if (tail != NEG_INF) total[relation].add(head + tail);
            }
        }
    zero = total[0].value();
    one = total[1].value();
}

void negate_costs(Problem& p, int j, vector<double>& table) {
    state_costs(p, j, table);
    for (double& value : table) value = -value;
}

}  // namespace

extern "C" {

// `(log Z_0(b), log Z_1(b))` for every boundary into `zero` and `one`.
int pbdp_posterior(int n, const int32_t* offsets, const int32_t* active, const double* c0,
                   const double* c1, int n_fragments, int interval, double* zero, double* one) {
    if (n < 2) return 0;
    Problem p{n, offsets, active, c0, c1, vector<int32_t>(size_t(n_fragments), -1)};
    Scratch s;
    const int step = interval > 0 ? interval : std::max(1, int(std::ceil(std::sqrt(double(n)))));
    const int segments = (n + step - 1) / step;
    vector<vector<double>> checkpoint(static_cast<size_t>(segments));
    vector<double> current, next;
    negate_costs(p, 0, current);
    checkpoint[0] = current;
    for (int j = 1; j < n; ++j) {
        sum_forward_step(p, s, current, j, next);
        current.swap(next);
        if (j % step == 0) checkpoint[size_t(j / step)] = current;
    }
    vector<double> backward, scratch;
    negate_costs(p, n - 1, backward);
    vector<vector<double>> table(static_cast<size_t>(step));
    for (int seg = segments - 1; seg >= 0; --seg) {
        const int a = seg * step;
        const int b = std::min(a + step - 1, n - 1);
        table[0] = checkpoint[size_t(seg)];
        for (int j = a + 1; j <= b; ++j) sum_forward_step(p, s, table[size_t(j - 1 - a)], j, table[size_t(j - a)]);
        for (int j = b; j >= a; --j) {
            if (j < n - 1) {
                sum_boundary(p, s, table[size_t(j - a)], backward, j, zero[j], one[j]);
                sum_backward_step(p, s, backward, j, scratch);
                backward.swap(scratch);
            }
        }
        checkpoint[size_t(seg)].clear();
        checkpoint[size_t(seg)].shrink_to_fit();
    }
    return 0;
}


// Returns 0 on success, 1 for a malformed matrix, 2 for a broken traceback.
// `zero` and `one` receive `boundary_costs` for every boundary when
// `want_margins` is nonzero; they may be null otherwise.
int pbdp_run(int n, const int32_t* offsets, const int32_t* active, const double* c0,
             const double* c1, int n_fragments, int want_margins, int interval,
             double* cost, int8_t* haplotype, int8_t* assignment, double* zero, double* one) {
    for (int f = 0; f < n_fragments; ++f) assignment[f] = 0;
    if (n == 0) {
        *cost = 0.0;
        return 0;
    }
    Problem p{n, offsets, active, c0, c1, vector<int32_t>(size_t(n_fragments), -1)};
    Scratch s;
    const int step = interval > 0 ? interval : std::max(1, int(std::ceil(std::sqrt(double(n)))));
    const int segments = (n + step - 1) / step;

    vector<vector<double>> checkpoint(static_cast<size_t>(segments));
    vector<double> current, next;
    state_costs(p, 0, current);
    checkpoint[0] = current;
    for (int j = 1; j < n; ++j) {
        forward_step(p, s, current, j, next);
        current.swap(next);
        if (j % step == 0) checkpoint[size_t(j / step)] = current;
    }

    size_t best_state = 0;
    for (size_t state = 1; state < current.size(); ++state)
        if (current[state] < current[best_state]) best_state = state;
    if (current[best_state] == INF) return 1;
    *cost = current[best_state];

    vector<int64_t> states(size_t(n), -1);
    states[size_t(n - 1)] = int64_t(best_state);
    vector<double> backward, scratch;
    if (want_margins) state_costs(p, n - 1, backward);

    vector<vector<double>> table(static_cast<size_t>(step));
    for (int seg = segments - 1; seg >= 0; --seg) {
        const int a = seg * step;
        const int b = std::min(a + step - 1, n - 1);
        table[0] = checkpoint[size_t(seg)];
        for (int j = a + 1; j <= b; ++j) forward_step(p, s, table[size_t(j - 1 - a)], j, table[size_t(j - a)]);
        if (b < n - 1) {
            states[size_t(b)] = parent_of(p, s, table[size_t(b - a)], b + 1, states[size_t(b + 1)]);
            if (states[size_t(b)] < 0) return 2;
        }
        for (int j = b; j >= a; --j) {
            if (want_margins && j < n - 1) {
                boundary(p, s, table[size_t(j - a)], backward, j, zero[j], one[j]);
                backward_step(p, s, backward, j, scratch);
                backward.swap(scratch);
            }
            if (j > a) {
                states[size_t(j - 1)] = parent_of(p, s, table[size_t(j - 1 - a)], j, states[size_t(j)]);
                if (states[size_t(j - 1)] < 0) return 2;
            }
        }
        checkpoint[size_t(seg)].clear();
        checkpoint[size_t(seg)].shrink_to_fit();
    }

    for (int j = 0; j < n; ++j) {
        haplotype[j] = int8_t(states[size_t(j)] & 1);
        const int64_t mask = states[size_t(j)] >> 1;
        for (int slot = 0; slot < p.width(j); ++slot)
            assignment[p.active[p.offsets[j] + slot]] = int8_t((mask >> slot) & 1);
    }
    return 0;
}

}  // extern "C"
