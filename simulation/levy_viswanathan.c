/*
levy_viswanathan.c — VISWANATHAN CONVENTION (flight interrupt + restart from target)

Rule (Viswanathan et al., 1999):
  The searcher flies along a straight trajectory of length L.
  If a target is detected within rv, the searcher:
    1) Moves to the target center (visit is counted).
    2) Flight TERMINATES immediately — remaining planned distance is discarded.
    3) A new flight begins from the target position with a freshly sampled
       direction theta and length L.

  This is the canonical "interrupted search" model from the Lévy foraging literature.

  Supports destructive (--remove-prob 1) and non-destructive (--remove-prob 0) modes.

Build:
  gcc -O3 -std=c11 -march=native -o levy_viswanathan levy_viswanathan.c -lm
*/

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ========================= Basic geometry =========================
typedef struct { double x, y; } Vec2;

static inline double dot2(Vec2 a, Vec2 b) { return a.x*b.x + a.y*b.y; }
static inline Vec2   add2(Vec2 a, Vec2 b) { Vec2 r={a.x+b.x,a.y+b.y}; return r; }
static inline Vec2   sub2(Vec2 a, Vec2 b) { Vec2 r={a.x-b.x,a.y-b.y}; return r; }
static inline Vec2   mul2(Vec2 a, double s){ Vec2 r={a.x*s,a.y*s}; return r; }
static inline double norm2(Vec2 a) { return a.x*a.x + a.y*a.y; }

static inline double dist(Vec2 a, Vec2 b) {
    double n2 = norm2(sub2(a,b));
    return (n2 >= 0.0) ? sqrt(n2) : 0.0;
}
static inline int i_floor_div(double x, double s) { return (int)floor(x / s); }

// ========================= RNG (xorshift64*) =========================
typedef struct { uint64_t s; } rng64_t;

static inline uint64_t splitmix64_next(uint64_t *x){
    uint64_t z = (*x += 0x9E3779B97f4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static inline uint64_t mix64(uint64_t x){ return splitmix64_next(&x); }

static inline void rng_seed(rng64_t* r, uint64_t seed){
    uint64_t x = seed;
    r->s = splitmix64_next(&x);
    if(!r->s) r->s = 0x9E3779B97f4A7C15ULL;
}
static inline uint64_t rng_next(rng64_t* r){
    uint64_t x=r->s;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    r->s = x;
    return x * 0x2545F4914F6CDD1DULL;
}
static inline double rng01(rng64_t* r){
    return (rng_next(r) >> 11) * (1.0/9007199254740992.0);
}
static inline double rng_angle(rng64_t* r){
    return rng01(r) * 2.0 * M_PI;
}

// ========================= Poisson sampler =========================
static int pois_knuth(rng64_t* r, double lam){
    double L = exp(-lam), p = 1.0;
    int k = 0;
    do { ++k; p *= rng01(r); } while(p > L);
    return k - 1;
}
static int pois_sample(rng64_t* r, double lam){
    if(lam <= 0.0) return 0;
    if(lam <= 30.0) return pois_knuth(r, lam);
    // normal approx for large lam
    double u1=fmax(1e-12,rng01(r)), u2=rng01(r);
    double z=sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
    int k=(int) llround(lam + sqrt(lam)*z);
    return (k < 0) ? 0 : k;
}

// ========================= Truncated power-law =========================
// pdf ∝ l^{-mu} on [lmin,lmax]
static double sample_powerlaw(rng64_t* r, double mu, double lmin, double lmax){
    double u = rng01(r);
    if(fabs(mu-1.0) < 1e-12) {
        return exp(log(lmin) + u*(log(lmax)-log(lmin)));
    }
    double a = 1.0 - mu;
    double lmin_a = pow(lmin,a), lmax_a = pow(lmax,a);
    return pow(lmin_a + u*(lmax_a-lmin_a), 1.0/a);
}

// ========================= Id set (hash set) =========================
// open addressing, 0 key means empty. Insert-only.
typedef struct {
    uint64_t* keys;
    size_t cap;   // power of 2
    size_t n;
} IdSet;

static void idset_init(IdSet* s, size_t cap_pow2){
    s->cap = cap_pow2;
    s->n = 0;
    s->keys = (uint64_t*)calloc(s->cap, sizeof(uint64_t));
    if(!s->keys){ fprintf(stderr,"OOM idset_init\n"); exit(1); }
}
static void idset_free(IdSet* s){
    free(s->keys);
    s->keys = NULL;
    s->cap = s->n = 0;
}
static inline size_t idset_probe(IdSet* s, uint64_t key){
    size_t mask = s->cap - 1;
    size_t i = (size_t)mix64(key) & mask;
    while(1){
        uint64_t cur = s->keys[i];
        if(cur == 0 || cur == key) return i;
        i = (i + 1) & mask;
    }
}
static int idset_contains(IdSet* s, uint64_t key){
    if(key == 0) return 0;
    size_t i = idset_probe(s, key);
    return s->keys[i] == key;
}
static void idset_rehash(IdSet* s, size_t new_cap){
    uint64_t* old = s->keys;
    size_t old_cap = s->cap;

    s->keys = (uint64_t*)calloc(new_cap, sizeof(uint64_t));
    if(!s->keys){ fprintf(stderr,"OOM idset_rehash\n"); exit(1); }
    s->cap = new_cap;
    s->n = 0;

    for(size_t i=0;i<old_cap;i++){
        uint64_t k = old[i];
        if(k){
            size_t j = idset_probe(s, k);
            s->keys[j] = k;
            s->n++;
        }
    }
    free(old);
}
static void idset_insert(IdSet* s, uint64_t key){
    if(key == 0) key = 1;
    if((s->n + 1) * 10 >= s->cap * 7){
        idset_rehash(s, s->cap * 2);
    }
    size_t i = idset_probe(s, key);
    if(s->keys[i] == 0){
        s->keys[i] = key;
        s->n++;
    }
}

// ========================= Chunked stationary PPP world =========================
typedef struct {
    double x, y;
    uint64_t id;
} Target;

typedef struct {
    int cx, cy;
    int n_targets;
    Target* targets;
    unsigned char* alive;
} Chunk;

static inline uint64_t chunk_key(int cx, int cy){
    return ((uint64_t)(uint32_t)cx << 32) | (uint64_t)(uint32_t)cy;
}

// state: 0 empty, 1 occupied, 2 tombstone
typedef struct {
    uint8_t state;
    uint64_t key;
    uint64_t last_used;
    Chunk chunk;
} ChunkSlot;

typedef struct {
    double chunk_size;
    double rho;
    uint64_t world_seed;
    uint64_t tick;

    ChunkSlot* slots;
    size_t table_cap;     // power of 2
    size_t n_active;
    size_t max_active;

    IdSet* removed;       // per-run removal state (destructive)
} ChunkManager;

static void chunk_free(Chunk* c){
    free(c->targets);
    free(c->alive);
    c->targets = NULL;
    c->alive = NULL;
    c->n_targets = 0;
}

static void cm_init(ChunkManager* cm,
                    double chunk_size, double rho,
                    uint64_t world_seed,
                    size_t table_cap_pow2,
                    size_t max_active,
                    IdSet* removed){
    cm->chunk_size = chunk_size;
    cm->rho = rho;
    cm->world_seed = world_seed;
    cm->tick = 1;

    cm->table_cap = table_cap_pow2;
    cm->slots = (ChunkSlot*)calloc(cm->table_cap, sizeof(ChunkSlot));
    if(!cm->slots){ fprintf(stderr,"OOM cm_init\n"); exit(1); }

    cm->n_active = 0;
    cm->max_active = max_active;
    cm->removed = removed;
}

static void cm_free(ChunkManager* cm){
    for(size_t i=0;i<cm->table_cap;i++){
        if(cm->slots[i].state == 1){
            chunk_free(&cm->slots[i].chunk);
        }
    }
    free(cm->slots);
    cm->slots = NULL;
    cm->table_cap = cm->n_active = cm->max_active = 0;
}

static inline size_t cm_find_slot(ChunkManager* cm, uint64_t key, int* found){
    size_t mask = cm->table_cap - 1;
    size_t i = (size_t)mix64(key) & mask;
    size_t first_tomb = (size_t)-1;

    while(1){
        uint8_t st = cm->slots[i].state;
        if(st == 0){
            *found = 0;
            return (first_tomb != (size_t)-1) ? first_tomb : i;
        }
        if(st == 1 && cm->slots[i].key == key){
            *found = 1;
            return i;
        }
        if(st == 2 && first_tomb == (size_t)-1){
            first_tomb = i;
        }
        i = (i + 1) & mask;
    }
}

static void cm_generate_chunk(ChunkManager* cm, Chunk* out, int cx, int cy){
    uint64_t key = chunk_key(cx, cy);
    uint64_t seed = mix64(cm->world_seed ^ key);

    rng64_t rng;
    rng_seed(&rng, seed);

    double area = cm->chunk_size * cm->chunk_size;
    int n = pois_sample(&rng, cm->rho * area);
    if(n < 0) n = 0;

    out->cx = cx;
    out->cy = cy;
    out->n_targets = n;

    out->targets = (n > 0) ? (Target*)malloc((size_t)n * sizeof(Target)) : NULL;
    out->alive   = (n > 0) ? (unsigned char*)malloc((size_t)n * sizeof(unsigned char)) : NULL;
    if((n>0) && (!out->targets || !out->alive)){
        fprintf(stderr,"OOM cm_generate_chunk n=%d\n", n);
        exit(1);
    }

    for(int i=0;i<n;i++){
        double ux = rng01(&rng);
        double uy = rng01(&rng);
        out->targets[i].x = (double)cx * cm->chunk_size + ux * cm->chunk_size;
        out->targets[i].y = (double)cy * cm->chunk_size + uy * cm->chunk_size;

        uint64_t id = mix64(key ^ (uint64_t)i * 0x9E3779B97f4A7C15ULL);
        if(id == 0) id = 1;
        out->targets[i].id = id;

        if(cm->removed){
            out->alive[i] = idset_contains(cm->removed, id) ? 0 : 1;
        } else {
            out->alive[i] = 1;
        }
    }
}

static void cm_evict_one(ChunkManager* cm){
    uint64_t best_tick = UINT64_MAX;
    size_t best_i = (size_t)-1;

    for(size_t i=0;i<cm->table_cap;i++){
        if(cm->slots[i].state == 1){
            if(cm->slots[i].last_used < best_tick){
                best_tick = cm->slots[i].last_used;
                best_i = i;
            }
        }
    }
    if(best_i == (size_t)-1) return;

    chunk_free(&cm->slots[best_i].chunk);
    cm->slots[best_i].state = 2;     // tombstone
    cm->slots[best_i].key = 0;
    cm->slots[best_i].last_used = 0;
    if(cm->n_active > 0) cm->n_active--;
}

static Chunk* cm_get_chunk(ChunkManager* cm, int cx, int cy){
    uint64_t key = chunk_key(cx, cy);
    int found = 0;
    size_t idx = cm_find_slot(cm, key, &found);

    if(found){
        cm->slots[idx].last_used = cm->tick++;
        return &cm->slots[idx].chunk;
    }

    if(cm->n_active >= cm->max_active){
        cm_evict_one(cm);
        idx = cm_find_slot(cm, key, &found);
        if(found){
            cm->slots[idx].last_used = cm->tick++;
            return &cm->slots[idx].chunk;
        }
    }

    cm->slots[idx].state = 1;
    cm->slots[idx].key = key;
    cm->slots[idx].last_used = cm->tick++;
    cm_generate_chunk(cm, &cm->slots[idx].chunk, cx, cy);
    cm->n_active++;
    return &cm->slots[idx].chunk;
}

// ========================= Parameters =========================
typedef struct {
    // environment
    double rho;          // target density (PPP, targets per unit area)
    double rv;           // detection radius
    double chunk_size;   // chunk side length

    // levy
    double mu_min, mu_max, mu_step;
    double lmin, lmax;   // truncation for step lengths

    // motion / detection discretization
    double dx;           // segment length for scanning

    // experiment size
    int flights;         // flights per run
    int reps;            // runs per mu
    int seed;            // base seed

    // destructive toggle
    double remove_prob;  // probability to remove visited target (0..1)

    // output
    const char* out_csv;

    // scoring (optional)
    double r_target;        // reward per unique capture
    double r_step;          // penalty per step
    double score_dx;        // step length for converting distance->steps
    int score_per_segment;  // if 1: steps=walks; else steps=distance/score_dx

    // convenience: expected targets per chunk
    double chunk_targ;      // if >0: rho = chunk_targ/(chunk_size^2)
} Params;

static inline double mean_free_path_2d(double rho, double rv){
    // Viswanathan lambda = 1/(2*rho*rv)
    if(rho <= 0.0 || rv <= 0.0) return NAN;
    return 1.0 / (2.0 * rv * rho);
}

// ========================= Run metrics (per rep) =========================
typedef struct {
    int captures_unique;     // unique target sites visited
    int visits;              // ALL visits (each arrival to a target counts)
    double total_distance;
    double total_time;       // speed=1

    int walks;               // number of scanned segments executed
    int num_flights;

    // first capture time (for MFPT-like diagnostics)
    int has_first;
    double first_capture_time;

    // inter-visit timing (optional)
    double last_visit_time;
    double sum_intervals;
    int cnt_intervals;
} RunMetrics;

static void init_metrics(RunMetrics* m){
    m->captures_unique = 0;
    m->visits = 0;
    m->total_distance = 0.0;
    m->total_time = 0.0;
    m->walks = 0;
    m->num_flights = 0;

    m->has_first = 0;
    m->first_capture_time = NAN;

    m->last_visit_time = NAN;
    m->sum_intervals = 0.0;
    m->cnt_intervals = 0;
}

static inline void move_to(RunMetrics* m, Vec2* pos, Vec2 new_pos){
    double d = dist(*pos, new_pos);
    if(!isfinite(d) || d < 0.0) d = 0.0;
    m->total_distance += d;
    m->total_time += d; // speed=1
    *pos = new_pos;
}

static inline void record_visit_time(RunMetrics* m, double t){
    if(!m->has_first){
        m->has_first = 1;
        m->first_capture_time = t;
        m->last_visit_time = t;
    } else {
        m->sum_intervals += (t - m->last_visit_time);
        m->cnt_intervals += 1;
        m->last_visit_time = t;
    }
}

// ========================= Detection helpers =========================
typedef struct {
    Chunk* chunk;
    int idx;
    double dist_to_target;
} TargetRef;

static int find_nearest_visible_target(ChunkManager* cm,
                                      Vec2 pos,
                                      double rv,
                                      uint64_t last_id,
                                      double ignore_last_within,
                                      TargetRef* out_ref){
    double rv2 = rv*rv;

    double minx = pos.x - rv;
    double maxx = pos.x + rv;
    double miny = pos.y - rv;
    double maxy = pos.y + rv;

    int cx0 = i_floor_div(minx, cm->chunk_size);
    int cx1 = i_floor_div(maxx, cm->chunk_size);
    int cy0 = i_floor_div(miny, cm->chunk_size);
    int cy1 = i_floor_div(maxy, cm->chunk_size);

    double best_d2 = INFINITY;
    Chunk* best_chunk = NULL;
    int best_idx = -1;

    for(int cy=cy0; cy<=cy1; cy++){
        for(int cx=cx0; cx<=cx1; cx++){
            Chunk* ch = cm_get_chunk(cm, cx, cy);
            for(int i=0;i<ch->n_targets;i++){
                if(!ch->alive[i]) continue;
                Target t = ch->targets[i];
                Vec2 d = sub2((Vec2){t.x,t.y}, pos);
                double d2 = d.x*d.x + d.y*d.y;

                if(d2 <= rv2){
                    double dlen = sqrt(d2);
                    if(t.id == last_id && dlen < ignore_last_within) continue;
                    if(d2 < best_d2){
                        best_d2 = d2;
                        best_chunk = ch;
                        best_idx = i;
                    }
                }
            }
        }
    }

    if(best_idx >= 0){
        out_ref->chunk = best_chunk;
        out_ref->idx = best_idx;
        out_ref->dist_to_target = sqrt(best_d2);
        return 1;
    }
    return 0;
}

// segment entry into circle of radius rv around c, starting at pos, direction dir (unit),
// for segment length seg_len. Returns earliest s in [0,seg_len].
static int segment_entry_time(Vec2 pos, Vec2 dir, double seg_len, Vec2 c, double rv, double* out_s){
    Vec2 rel = sub2(pos, c);
    double b = dot2(rel, dir);
    double c0 = dot2(rel, rel) - rv*rv;
    double disc = b*b - c0;
    if(disc < 0.0) return 0;
    double sqrt_disc = sqrt(disc);
    double s1 = -b - sqrt_disc;
    double s2 = -b + sqrt_disc;

    if(s2 < 0.0) return 0;
    double s = (s1 < 0.0) ? 0.0 : s1;
    if(s > seg_len) return 0;
    *out_s = s;
    return 1;
}

// Find earliest detection along the segment, over any target disks.
static int find_first_detection_in_segment(ChunkManager* cm,
                                          Vec2 pos,
                                          Vec2 dir, double seg_len,
                                          double rv,
                                          uint64_t last_id,
                                          double ignore_last_within,
                                          double* out_s_hit,
                                          Vec2* out_det_point){
    double minx = fmin(pos.x, pos.x + dir.x*seg_len) - rv;
    double maxx = fmax(pos.x, pos.x + dir.x*seg_len) + rv;
    double miny = fmin(pos.y, pos.y + dir.y*seg_len) - rv;
    double maxy = fmax(pos.y, pos.y + dir.y*seg_len) + rv;

    int cx0 = i_floor_div(minx, cm->chunk_size);
    int cx1 = i_floor_div(maxx, cm->chunk_size);
    int cy0 = i_floor_div(miny, cm->chunk_size);
    int cy1 = i_floor_div(maxy, cm->chunk_size);

    double best_s = INFINITY;

    for(int cy=cy0; cy<=cy1; cy++){
        for(int cx=cx0; cx<=cx1; cx++){
            Chunk* ch = cm_get_chunk(cm, cx, cy);
            for(int i=0;i<ch->n_targets;i++){
                if(!ch->alive[i]) continue;
                Target t = ch->targets[i];

                if(t.id == last_id){
                    double d0 = dist(pos, (Vec2){t.x,t.y});
                    if(d0 < ignore_last_within) continue;
                }

                double s_hit;
                if(segment_entry_time(pos, dir, seg_len, (Vec2){t.x,t.y}, rv, &s_hit)){
                    if(s_hit < best_s){
                        best_s = s_hit;
                    }
                }
            }
        }
    }

    if(!isfinite(best_s)) return 0;
    *out_s_hit = best_s;
    *out_det_point = add2(pos, mul2(dir, best_s));
    return 1;
}

// ========================= Destructive removal =========================
static void remove_target_if_needed(rng64_t* rng, const Params* P, ChunkManager* cm, Chunk* ch, int idx){
    if(P->remove_prob <= 0.0) return;
    if(rng01(rng) < P->remove_prob){
        uint64_t id = ch->targets[idx].id;
        ch->alive[idx] = 0;
        if(cm->removed) {
            idset_insert(cm->removed, id);
        }
    }
}

// ========================= One run (correct flight termination) =========================
static void run_one_series(uint64_t series_seed, const Params* P, double mu, RunMetrics* metrics){
    rng64_t rng;
    rng_seed(&rng, series_seed);

    // removed-set only for destructive mode
    IdSet removed;
    IdSet* removed_ptr = NULL;
    if(P->remove_prob > 0.0){
        idset_init(&removed, 1u<<18);
        removed_ptr = &removed;
    }

    // visited-set for UNIQUE target sites
    IdSet visited;
    idset_init(&visited, 1u<<18);

    uint64_t world_seed = mix64(series_seed ^ 0xD6E8FEB86659FD93ULL);

    ChunkManager cm;
    cm_init(&cm, P->chunk_size, P->rho, world_seed,
            1u<<15,   // table cap = 32768
            4096,     // active chunk cache
            removed_ptr);

    Vec2 pos = (Vec2){0.0, 0.0};
    uint64_t last_id = 0;

    // After visiting a target, ignore it for the ENTIRE next flight.
    // This prevents "ping-pong" trapping where the agent immediately
    // re-detects the same target after a short flight.
    // Viswanathan's rule: "chooses a new direction and a new distance l_{j+1}"
    // implies the forager needs one full flight before re-detecting.
    int just_visited = 0;

    for(int k=0; k<P->flights; ++k){
        metrics->num_flights++;

        double theta = rng_angle(&rng);
        Vec2 dir = (Vec2){cos(theta), sin(theta)};
        double planned = sample_powerlaw(&rng, mu, P->lmin, P->lmax);

        double flight_left = planned;
        int captured = 0;

        // If we just visited last_id, ignore it for this ENTIRE flight
        double ignore_last_within = just_visited ? 1e30 : 0.0;
        just_visited = 0;  // reset for next iteration

        while(flight_left > 0.0){
            double seg_len = fmin(P->dx, flight_left);

            metrics->walks++;

            double s_hit;
            Vec2 det_point;

            if(find_first_detection_in_segment(&cm, pos, dir, seg_len, P->rv,
                                               last_id, ignore_last_within,
                                               &s_hit, &det_point)) {
                // Move to detection point along the flight
                move_to(metrics, &pos, det_point);

                // From detection point, select nearest visible target and go to its center
                TargetRef ref;
                if(!find_nearest_visible_target(&cm, pos, P->rv, last_id, ignore_last_within, &ref)){
                    // Rare numerical/geometry mismatch: if cannot resolve target, just move to end of this segment
                    Vec2 endp = add2(pos, mul2(dir, seg_len - s_hit));
                    move_to(metrics, &pos, endp);
                    flight_left -= seg_len;
                    continue;
                }

                Target tgt = ref.chunk->targets[ref.idx];
                Vec2 tgt_pos = (Vec2){tgt.x, tgt.y};

                // Move to target center (this is part of the path/time)
                move_to(metrics, &pos, tgt_pos);

                // Count VISIT (every arrival counts)
                metrics->visits++;
                record_visit_time(metrics, metrics->total_time);

                // Count UNIQUE site
                if(!idset_contains(&visited, tgt.id)){
                    idset_insert(&visited, tgt.id);
                    metrics->captures_unique++;
                }

                // Potentially remove target (destructive)
                remove_target_if_needed(&rng, P, &cm, ref.chunk, ref.idx);

                last_id = tgt.id;
                just_visited = 1;  // ignore last_id for the ENTIRE next flight

                // Viswanathan rule: terminate flight, start new one from target
                break;
            } else {
                // No detection in this segment: move to segment end
                Vec2 endp = add2(pos, mul2(dir, seg_len));
                move_to(metrics, &pos, endp);
                flight_left -= seg_len;
            }
        }
        // Next k => new (theta, planned)
    }

    cm_free(&cm);
    idset_free(&visited);
    if(removed_ptr) idset_free(removed_ptr);
}

// ========================= Stats helpers =========================
typedef struct {
    double mean, std, var;
    double p10, p90, min, max;
} Stats7;

static int dbl_cmp(const void* a, const void* b){
    double da=*(const double*)a, db=*(const double*)b;
    return (da<db)?-1:((da>db)?1:0);
}
static double percentile(double* arr, int n, double p){
    if(n<=0) return NAN;
    if(p<=0) return arr[0];
    if(p>=1) return arr[n-1];
    double pos = p*(n-1);
    int i=(int)floor(pos);
    double frac=pos - i;
    if(i>=n-1) return arr[n-1];
    return arr[i]*(1.0-frac) + arr[i+1]*frac;
}
static int compact_finite(double* dst, const double* src, int n){
    int m = 0;
    for(int i=0;i<n;i++){
        if(isfinite(src[i])) dst[m++] = src[i];
    }
    return m;
}
static void stats7_from_array(double* arr, int n, Stats7* out){
    if(n <= 0){
        out->mean = out->std = out->var = NAN;
        out->p10 = out->p90 = out->min = out->max = NAN;
        return;
    }
    double sum = 0.0, sum2 = 0.0;
    for(int i=0;i<n;i++){
        sum += arr[i];
        sum2 += arr[i]*arr[i];
    }
    out->mean = sum / (double)n;
    out->var  = fmax(0.0, sum2/(double)n - out->mean*out->mean);
    out->std  = sqrt(out->var);

    qsort(arr, (size_t)n, sizeof(double), dbl_cmp);
    out->min = arr[0];
    out->max = arr[n-1];
    out->p10 = percentile(arr, n, 0.1);
    out->p90 = percentile(arr, n, 0.9);
}

// ========================= Aggregation =========================
typedef struct {
    double rho, rv, mu, lmin, lmax, dx;
    double chunk_size, remove_prob;
    int flights, seed, reps;

    Stats7 unique_caps;  // captures_unique per rep
    Stats7 visits;       // visits per rep

    double mean_distance, mean_time;

    double mfpt;               // mean first capture time (if none => total_time)
    double mean_time_between;  // mean inter-visit time across intervals

    // Efficiency:
    // eta_unique = unique_caps/distance, eta_visits = visits/distance
    double eta_unique_global;
    double eta_visits_global;

    // Normalized rate:
    // normalized_rate_visits_global = lambda * eta_visits_global   (KEY METRIC)
    double normalized_rate_unique_global;
    double normalized_rate_visits_global;

    // Per-rep distributions (optional)
    Stats7 eta_unique;
    Stats7 eta_visits;
    Stats7 normalized_rate_unique;
    Stats7 normalized_rate_visits;

    // Score (optional)
    Stats7 score;
} AggregatedMetrics;

static inline double compute_score(const Params* P, const RunMetrics* m){
    double steps;
    if(P->score_per_segment){
        steps = (double)m->walks;
    } else {
        steps = (P->score_dx > 0.0) ? (m->total_distance / P->score_dx) : NAN;
    }
    return P->r_target * (double)m->captures_unique + P->r_step * steps;
}

static void aggregate_metrics(const RunMetrics* runs, int n_reps, const Params* P, double mu, AggregatedMetrics* out){
    out->rho = P->rho; out->rv = P->rv; out->mu = mu;
    out->lmin = P->lmin; out->lmax = P->lmax; out->dx = P->dx;
    out->chunk_size = P->chunk_size;
    out->remove_prob = P->remove_prob;
    out->flights = P->flights;
    out->seed = P->seed;
    out->reps = n_reps;

    double sum_dist=0.0, sum_time=0.0;

    double sum_mfpt=0.0; int cnt_mfpt=0;
    double sum_intervals=0.0; int cnt_intervals=0;

    double tot_unique=0.0, tot_visits=0.0, tot_dist=0.0;

    double lambda = mean_free_path_2d(P->rho, P->rv);

    double* uniq_raw  = (double*)malloc((size_t)n_reps * sizeof(double));
    double* vis_raw   = (double*)malloc((size_t)n_reps * sizeof(double));
    double* eta_u_raw = (double*)malloc((size_t)n_reps * sizeof(double));
    double* eta_v_raw = (double*)malloc((size_t)n_reps * sizeof(double));
    double* nr_u_raw  = (double*)malloc((size_t)n_reps * sizeof(double));
    double* nr_v_raw  = (double*)malloc((size_t)n_reps * sizeof(double));
    double* score_raw = (double*)malloc((size_t)n_reps * sizeof(double));

    if(!uniq_raw || !vis_raw || !eta_u_raw || !eta_v_raw || !nr_u_raw || !nr_v_raw || !score_raw){
        fprintf(stderr,"OOM aggregate arrays\n"); exit(1);
    }

    for(int i=0;i<n_reps;i++){
        double uniq = (double)runs[i].captures_unique;
        double vis  = (double)runs[i].visits;
        double dist_i = runs[i].total_distance;
        double time_i = runs[i].total_time;

        uniq_raw[i] = uniq;
        vis_raw[i]  = vis;

        double eta_u = (dist_i > 0.0) ? (uniq / dist_i) : NAN;
        double eta_v = (dist_i > 0.0) ? (vis  / dist_i) : NAN;

        eta_u_raw[i] = eta_u;
        eta_v_raw[i] = eta_v;

        nr_u_raw[i] = (isfinite(lambda) && isfinite(eta_u)) ? (lambda * eta_u) : NAN;
        nr_v_raw[i] = (isfinite(lambda) && isfinite(eta_v)) ? (lambda * eta_v) : NAN;

        score_raw[i] = compute_score(P, &runs[i]);

        sum_dist += dist_i;
        sum_time += time_i;

        tot_unique += uniq;
        tot_visits += vis;
        tot_dist += dist_i;

        double mfpt_i = runs[i].has_first ? runs[i].first_capture_time : time_i;
        if(isfinite(mfpt_i)){ sum_mfpt += mfpt_i; cnt_mfpt++; }

        if(runs[i].cnt_intervals > 0){
            sum_intervals += runs[i].sum_intervals;
            cnt_intervals += runs[i].cnt_intervals;
        }
    }

    // stats over arrays (finite only)
    double* tmp = (double*)malloc((size_t)n_reps * sizeof(double));
    if(!tmp){ fprintf(stderr,"OOM tmp\n"); exit(1); }

    int m1 = compact_finite(tmp, uniq_raw, n_reps);
    stats7_from_array(tmp, m1, &out->unique_caps);

    int m2 = compact_finite(tmp, vis_raw, n_reps);
    stats7_from_array(tmp, m2, &out->visits);

    int m3 = compact_finite(tmp, eta_u_raw, n_reps);
    stats7_from_array(tmp, m3, &out->eta_unique);

    int m4 = compact_finite(tmp, eta_v_raw, n_reps);
    stats7_from_array(tmp, m4, &out->eta_visits);

    int m5 = compact_finite(tmp, nr_u_raw, n_reps);
    stats7_from_array(tmp, m5, &out->normalized_rate_unique);

    int m6 = compact_finite(tmp, nr_v_raw, n_reps);
    stats7_from_array(tmp, m6, &out->normalized_rate_visits);

    int m7 = compact_finite(tmp, score_raw, n_reps);
    stats7_from_array(tmp, m7, &out->score);

    free(tmp);

    out->mean_distance = sum_dist / (double)n_reps;
    out->mean_time     = sum_time / (double)n_reps;

    out->mfpt = (cnt_mfpt > 0) ? (sum_mfpt / (double)cnt_mfpt) : NAN;
    out->mean_time_between = (cnt_intervals > 0) ? (sum_intervals / (double)cnt_intervals) : NAN;

    out->eta_unique_global = (tot_dist > 0.0) ? (tot_unique / tot_dist) : NAN;
    out->eta_visits_global = (tot_dist > 0.0) ? (tot_visits / tot_dist) : NAN;

    out->normalized_rate_unique_global =
        (isfinite(lambda) && isfinite(out->eta_unique_global)) ? (lambda * out->eta_unique_global) : NAN;

    // KEY METRIC: this is what you should compare to Fig.2-style “l*h” curves (with visits, not unique).
    out->normalized_rate_visits_global =
        (isfinite(lambda) && isfinite(out->eta_visits_global)) ? (lambda * out->eta_visits_global) : NAN;

    free(uniq_raw);
    free(vis_raw);
    free(eta_u_raw);
    free(eta_v_raw);
    free(nr_u_raw);
    free(nr_v_raw);
    free(score_raw);
}

// ========================= CLI =========================
static int has_flag(int argc, char** argv, const char* name){
    for(int i=1;i<argc;i++) if(strcmp(argv[i],name)==0) return 1;
    return 0;
}
static double get_argd(int argc, char** argv, const char* name, double def){
    for(int i=1;i<argc-1;i++) if(strcmp(argv[i],name)==0) return atof(argv[i+1]);
    return def;
}
static const char* get_args(int argc, char** argv, const char* name, const char* def){
    for(int i=1;i<argc-1;i++) if(strcmp(argv[i],name)==0) return argv[i+1];
    return def;
}

static void usage(void){
    printf(
"Usage: ./levy_batch_infinite [options]\n"
"\n"
"Environment:\n"
"  --rho <float>          target density (2D PPP), targets per unit area (default 1e-3)\n"
"  --rv <float>           detection radius (default 1)\n"
"  --chunk-size <float>   chunk side length (default 500)\n"
"  --chunk-targ <float>   expected targets per chunk; sets rho=chunk_targ/chunk_size^2\n"
"\n"
"Lévy flight:\n"
"  --mu-min <f>           start mu (default 1.0)\n"
"  --mu-max <f>           end mu (default 3.0)\n"
"  --mu-step <f>          step for mu (default 0.1)\n"
"  --lmin <f>             min flight length (default auto = rv)\n"
"  --lmax <f>             max flight length (default auto = lambda=1/(2*rho*rv))\n"
"\n"
"Motion / detection:\n"
"  --dx <f>               segment length for scanning (default 50)\n"
"\n"
"Experiment size:\n"
"  --flights <int>        number of flights per run (default 250)\n"
"  --reps <int>           runs per mu (default 100)\n"
"  --seed <int>           base seed (default 42)\n"
"  --random-seed          use time-based seed\n"
"\n"
"Destructive mode:\n"
"  --remove-prob <f>      delete visited target prob [0..1] (default 0.0)\n"
"\n"
"Scoring (optional):\n"
"  --r-target <f>         reward per unique target (default 1.0)\n"
"  --r-step <f>           penalty per step (default -1e-3)\n"
"  --score-dx <f>         step length used for steps ~= distance/score_dx (default = dx)\n"
"  --score-per-segment    use steps = walks (segments) instead of distance/score_dx\n"
"\n"
"Output:\n"
"  --out <path>           output CSV (default results_infinite.csv)\n"
"  --append               append to CSV instead of overwrite\n"
    );
}

int main(int argc, char** argv){
    if(has_flag(argc,argv,"--help") || has_flag(argc,argv,"-h")){
        usage(); return 0;
    }

    Params P = {
        .rho = 1e-3,
        .rv  = 1.0,
        .chunk_size = 500.0,
        .chunk_targ = 0.0,

        .mu_min = 1.0,
        .mu_max = 3.0,
        .mu_step = 0.1,

        .lmin = 0.0,     // auto
        .lmax = 0.0,     // auto

        .dx = 50.0,

        .flights = 250,
        .reps = 100,
        .seed = 42,

        .remove_prob = 0.0,

        .out_csv = "results_infinite.csv",

        .r_target = 1.0,
        .r_step = -1e-3,
        .score_dx = 0.0,          // auto = dx
        .score_per_segment = 0
    };

    // parse
    P.rho = get_argd(argc,argv,"--rho",P.rho);
    P.rv  = get_argd(argc,argv,"--rv",P.rv);
    P.chunk_size = get_argd(argc,argv,"--chunk-size",P.chunk_size);
    P.chunk_targ = get_argd(argc,argv,"--chunk-targ",P.chunk_targ);

    P.mu_min = get_argd(argc,argv,"--mu-min",P.mu_min);
    P.mu_max = get_argd(argc,argv,"--mu-max",P.mu_max);
    P.mu_step = get_argd(argc,argv,"--mu-step",P.mu_step);

    P.lmin = get_argd(argc,argv,"--lmin",P.lmin);
    P.lmax = get_argd(argc,argv,"--lmax",P.lmax);

    P.dx   = get_argd(argc,argv,"--dx",P.dx);

    P.flights = (int)llround(get_argd(argc,argv,"--flights",(double)P.flights));
    P.reps    = (int)llround(get_argd(argc,argv,"--reps",(double)P.reps));
    P.seed    = (int)llround(get_argd(argc,argv,"--seed",(double)P.seed));

    P.remove_prob = get_argd(argc,argv,"--remove-prob",P.remove_prob);

    P.r_target = get_argd(argc,argv,"--r-target",P.r_target);
    P.r_step   = get_argd(argc,argv,"--r-step",P.r_step);
    P.score_dx = get_argd(argc,argv,"--score-dx",P.score_dx);

    P.out_csv = get_args(argc,argv,"--out",P.out_csv);

    if(has_flag(argc,argv,"--random-seed")){
        P.seed = (int)time(NULL) ^ (int)clock();
    }
    if(has_flag(argc,argv,"--score-per-segment")){
        P.score_per_segment = 1;
    }

    // validate
    if(P.chunk_size <= 0){ fprintf(stderr,"Bad chunk_size\n"); return 1; }
    if(P.chunk_targ > 0.0){
        P.rho = P.chunk_targ / (P.chunk_size * P.chunk_size);
    }
    if(P.rv <= 0 || P.rho < 0){ fprintf(stderr,"Bad rv/rho\n"); return 1; }
    if(P.dx <= 0){ fprintf(stderr,"Bad dx\n"); return 1; }
    if(P.mu_step <= 0 || P.mu_min > P.mu_max){ fprintf(stderr,"Bad mu range\n"); return 1; }
    if(P.flights <= 0 || P.reps <= 0){ fprintf(stderr,"Bad flights/reps\n"); return 1; }
    if(P.remove_prob < 0.0 || P.remove_prob > 1.0){ fprintf(stderr,"Bad remove_prob\n"); return 1; }

    if(P.score_dx <= 0.0) P.score_dx = P.dx;

    // Viswanathan defaults: lmin=rv, lmax=lambda
    if(P.lmin <= 0.0) P.lmin = P.rv;
    if(P.lmax <= 0.0){
        double lambda = mean_free_path_2d(P.rho, P.rv);
        if(!isfinite(lambda) || lambda <= P.lmin){
            fprintf(stderr,"Auto lmax failed: lambda=%g (rho=%g rv=%g)\n", lambda, P.rho, P.rv);
            return 1;
        }
        P.lmax = lambda;
    }
    if(!(P.lmax > P.lmin && P.lmin > 0)){
        fprintf(stderr,"Require lmax > lmin > 0 (got lmin=%g lmax=%g)\n", P.lmin, P.lmax);
        return 1;
    }

    int append_mode = has_flag(argc,argv,"--append");
    FILE* f = fopen(P.out_csv, append_mode ? "a" : "w");
    if(!f){ perror("fopen"); return 1; }

    if(!append_mode){
        fprintf(f,
            "rho,rv,mu,lmin,lmax,dx,chunk_size,remove_prob,flights,seed,reps,lambda,"
            "unique_caps_mean,unique_caps_std,unique_caps_var,unique_caps_p10,unique_caps_p90,unique_caps_min,unique_caps_max,"
            "visits_mean,visits_std,visits_var,visits_p10,visits_p90,visits_min,visits_max,"
            "mean_distance,mean_time,mfpt,mean_time_between,"
            "eta_unique_global,eta_visits_global,"
            "normalized_rate_unique_global,normalized_rate_visits_global,"
            "eta_unique_mean,eta_unique_std,eta_unique_min,eta_unique_max,"
            "eta_visits_mean,eta_visits_std,eta_visits_min,eta_visits_max,"
            "normalized_rate_unique_mean,normalized_rate_unique_std,normalized_rate_unique_min,normalized_rate_unique_max,"
            "normalized_rate_visits_mean,normalized_rate_visits_std,normalized_rate_visits_min,normalized_rate_visits_max,"
            "score_mean,score_std,score_min,score_max,"
            "r_target,r_step,score_dx,score_mode\n"
        );
    }

    double lambda = mean_free_path_2d(P.rho, P.rv);

    // robust mu-loop: integer index to avoid drift
    int ai = 0;
    for(;; ai++){
        double mu = P.mu_min + (double)ai * P.mu_step;
        if(mu > P.mu_max + 1e-12) break;

        RunMetrics* runs = (RunMetrics*)malloc((size_t)P.reps * sizeof(RunMetrics));
        if(!runs){ fprintf(stderr,"OOM runs\n"); fclose(f); return 1; }

        uint64_t rho_hash = (uint64_t)llround(P.rho * 1e15);
        uint64_t rv_hash  = (uint64_t)llround(P.rv  * 1e15);

        for(int r=0;r<P.reps;r++){
            init_metrics(&runs[r]);

            uint64_t s = (uint64_t)P.seed
                       ^ (uint64_t)ai * 0x9E3779B97f4A7C15ULL
                       ^ (uint64_t)r  * 0xBF58476D1CE4E5B9ULL
                       ^ rho_hash * 0x517CC1B727220A95ULL
                       ^ rv_hash  * 0x85EBCA6B;

            run_one_series(s, &P, mu, &runs[r]);
        }

        AggregatedMetrics agg;
        aggregate_metrics(runs, P.reps, &P, mu, &agg);

        const char* score_mode = P.score_per_segment ? "segments" : "distance";

        fprintf(f,
            "%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%d,%d,%d,%.12g,"
            "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
            "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
            "%.6f,%.6f,%.6f,%.6f,"
            "%.12g,%.12g,"
            "%.12g,%.12g,"
            "%.12g,%.12g,%.12g,%.12g,"
            "%.12g,%.12g,%.12g,%.12g,"
            "%.12g,%.12g,%.12g,%.12g,"
            "%.12g,%.12g,%.12g,%.12g,"
            "%.6f,%.6f,%.6f,%.6f,"
            "%.6g,%.6g,%.6g,%s\n",
            P.rho, P.rv, mu, P.lmin, P.lmax, P.dx, P.chunk_size, P.remove_prob, P.flights, P.seed, P.reps, lambda,

            agg.unique_caps.mean, agg.unique_caps.std, agg.unique_caps.var, agg.unique_caps.p10, agg.unique_caps.p90, agg.unique_caps.min, agg.unique_caps.max,

            agg.visits.mean, agg.visits.std, agg.visits.var, agg.visits.p10, agg.visits.p90, agg.visits.min, agg.visits.max,

            agg.mean_distance, agg.mean_time, agg.mfpt, agg.mean_time_between,

            agg.eta_unique_global, agg.eta_visits_global,

            agg.normalized_rate_unique_global, agg.normalized_rate_visits_global,

            agg.eta_unique.mean, agg.eta_unique.std, agg.eta_unique.min, agg.eta_unique.max,
            agg.eta_visits.mean, agg.eta_visits.std, agg.eta_visits.min, agg.eta_visits.max,

            agg.normalized_rate_unique.mean, agg.normalized_rate_unique.std, agg.normalized_rate_unique.min, agg.normalized_rate_unique.max,
            agg.normalized_rate_visits.mean, agg.normalized_rate_visits.std, agg.normalized_rate_visits.min, agg.normalized_rate_visits.max,

            agg.score.mean, agg.score.std, agg.score.min, agg.score.max,

            P.r_target, P.r_step, P.score_dx, score_mode
        );

        fprintf(stderr,
            "mu=%.3f  visits=%.2f±%.2f  eta_visits=%.6g±%.2g  (l*eta_visits)=%.6g±%.2g  (rho=%g rv=%g lmax=%g reps=%d)\n",
            mu,
            agg.visits.mean, agg.visits.std,
            agg.eta_visits.mean, agg.eta_visits.std,
            agg.normalized_rate_visits.mean, agg.normalized_rate_visits.std,
            P.rho, P.rv, P.lmax, P.reps
        );

        free(runs);
    }

    fclose(f);
    fprintf(stderr,"Saved CSV to %s\n", P.out_csv);
    return 0;
}
