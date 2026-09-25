/* Protocol tests: ballots (completeness, serialization, negative cases), share checks,
 * interpolation, aggregation and tally. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ballot.h"
#include "codec.h"
#include "shamir.h"
#include "agg.h"
#include "tally.h"
#include "evolve.h"

static int fails = 0;
#define CHECK(cond, ...) do { if (!(cond)) { fails++; printf("FAIL: " __VA_ARGS__); printf("\n"); } else if (verbose) { printf("  ok: " __VA_ARGS__); printf("\n"); } } while (0)
static int verbose = 1;

static void seed_from(uint8_t s[32], unsigned x) { memset(s, 0, 32); for (int i = 0; i < 4; i++) s[i] = (uint8_t)(x >> (8 * i)); }

static void test_ballots(int n, int t, int L, int w, int d, uint64_t q) {
    printf("ballots: n=%d t=%d L=%d w=%d d=%d logq=%d\n", n, t, L, w, d, 63 - __builtin_clzll(q) + 1);
    tv_params prm;
    CHECK(tv_params_init(&prm, n, t, L, w, d, q) == 0, "params");
    uint8_t s[32];
    seed_from(s, 1);
    tv_pub pub;
    CHECK(tv_pub_init(&pub, &prm, s) == 0, "public parameters");
    tv_akey *keys = calloc((size_t)n, sizeof(tv_akey));
    for (int k = 0; k < n; k++) seed_from(keys[k].key, 100 + (unsigned)k);

    /* parity-check matrix annihilates every sharing */
    {
        prg g;
        prg_init(&g, 9, s, 32);
        poly *vv = calloc((size_t)L, sizeof(poly)), *m = calloc((size_t)(n + 1) * L, sizeof(poly));
        for (int a = 0; a < L; a++) sample_uniform_poly(&vv[a], &g);
        shamir_share(m, vv, &prm, &g);
        int ok = 1;
        for (int l = 0; l <= n - t; l++)
            for (int a = 0; a < L; a++) {
                poly acc, tmp;
                poly_zero(&acc);
                for (int j = 0; j <= n; j++) { poly_scale(&tmp, &m[j * L + a], pub.H_m[l * (n + 1) + j]); poly_add(&acc, &acc, &tmp); }
                for (int i = 0; i < TV_N; i++) ok &= acc.c[i] == 0;
            }
        CHECK(ok, "H * shares = 0");
        /* any t shares interpolate to the secret */
        int T[16];
        for (int i = 0; i < t; i++) T[i] = n - i;
        uint64_t lam[16];
        shamir_lagrange(lam, T, t);
        for (int a = 0; a < L; a++) {
            poly acc, tmp;
            poly_zero(&acc);
            for (int i = 0; i < t; i++) { poly_scale(&tmp, &m[T[i] * L + a], to_mont(lam[i])); poly_add(&acc, &acc, &tmp); }
            ok &= memcmp(&acc, &vv[a], sizeof(poly)) == 0;
        }
        CHECK(ok, "Lagrange interpolation from %d shares", t);
        free(vv); free(m);
    }

    int votes[3][16];
    int nvotes = 0;
    /* admissible votes */
    if (w == TV_W_FREE) {
        for (int a = 0; a < L; a++) { votes[0][a] = 0; votes[1][a] = 1; votes[2][a] = a % 2; }
        nvotes = 3;
    } else {
        for (int x = 0; x < 2; x++) { for (int a = 0; a < L; a++) votes[x][a] = (a >= x && a < x + w) ? 1 : 0; }
        nvotes = 2;
    }
    uint8_t *buf = malloc(tv_ballot_maxbytes(&prm));
    for (int x = 0; x < nvotes; x++) {
        tv_ballot b, b2;
        tv_voter_secret sec;
        tv_prove_stats st;
        tv_ballot_alloc(&b, &prm);
        tv_ballot_alloc(&b2, &prm);
        tv_secret_alloc(&sec, &prm);
        seed_from(s, 1000 + (unsigned)x);
        CHECK(tv_vote(&b, &sec, &st, &pub, keys, 77 + (uint64_t)x, votes[x], s) == 0, "vote %d produced (%d attempts)", x, st.attempts);
        CHECK(tv_verify_ballot(&pub, &b) == 1, "ballot %d verifies", x);
        tv_sizes sz;
        size_t len = tv_ballot_encode(buf, tv_ballot_maxbytes(&prm), &pub, &b, &sz);
        CHECK(len > 0 && tv_ballot_decode(&b2, &pub, buf, len) == 0, "encode/decode (%zu bytes: com %zu, ct %zu, proof %zu)", len, sz.commitments, sz.ciphertexts, sz.proof);
        CHECK(tv_verify_ballot(&pub, &b2) == 1, "decoded ballot verifies");
        CHECK(tv_ballot_decode(&b2, &pub, buf, len - 10) != 0, "truncated encoding rejected");
        buf[len] = 0;
        CHECK(tv_ballot_decode(&b2, &pub, buf, len + 1) != 0, "encoding with an extra byte rejected");
        /* shares */
        int ok = 1;
        poly *mk = calloc((size_t)L, sizeof(poly));
        int64_t *rk = malloc(sizeof(int64_t) * (size_t)prm.mu * TV_N);
        for (int k = 1; k <= n; k++) {
            ok &= tv_check_share(&pub, &keys[k - 1], k, &b, mk, rk) == 1;
            for (int a = 0; a < L; a++) ok &= memcmp(&mk[a], &sec.m[k * L + a], sizeof(poly)) == 0;
        }
        CHECK(ok, "every authority opens its share");
        CHECK(tv_check_share(&pub, &keys[0], 2, &b, mk, rk) == 0, "wrong key -> invalid opening (complaint)");
        /* negative cases: each tampering must be rejected */
        tv_ballot_decode(&b2, &pub, buf, len); b2.z[5] += 1;
        CHECK(!tv_verify_ballot(&pub, &b2), "reject: modified response z");
        tv_ballot_decode(&b2, &pub, buf, len); b2.or_r[7] -= 1;
        CHECK(!tv_verify_ballot(&pub, &b2), "reject: modified OR response");
        tv_ballot_decode(&b2, &pub, buf, len); b2.chat[3] ^= 1;
        CHECK(!tv_verify_ballot(&pub, &b2), "reject: modified challenge");
        tv_ballot_decode(&b2, &pub, buf, len); b2.c[prm.d].c[0] = mod_add(b2.c[prm.d].c[0], 1);
        CHECK(!tv_verify_ballot(&pub, &b2), "reject: commitment to vote + 1");
        tv_ballot_decode(&b2, &pub, buf, len); b2.c[prm.rows + prm.d].c[0] = mod_add(b2.c[prm.rows + prm.d].c[0], 1);
        CHECK(!tv_verify_ballot(&pub, &b2), "reject: inconsistent share");
        tv_ballot_decode(&b2, &pub, buf, len); b2.e[0] ^= 1;
        CHECK(!tv_verify_ballot(&pub, &b2), "reject: modified ciphertext");
        tv_ballot_decode(&b2, &pub, buf, len); b2.id ^= 1;
        CHECK(!tv_verify_ballot(&pub, &b2), "reject: other voter identity");
        tv_ballot_decode(&b2, &pub, buf, len);
        { challenge f = b2.f0[0]; f.sgn[0] = (int8_t)-f.sgn[0]; b2.f0[0] = f; }
        CHECK(!tv_verify_ballot(&pub, &b2), "reject: modified OR challenge f_0");
        tv_ballot_decode(&b2, &pub, buf, len);
        for (size_t i = 0; i < (size_t)(n + 1) * prm.mu * TV_N; i++) b2.z[i] *= 3;
        CHECK(!tv_verify_ballot(&pub, &b2), "reject: responses above the norm bound");
        free(mk); free(rk);
        tv_ballot_free(&b); tv_ballot_free(&b2); tv_secret_free(&sec);
    }
    /* inadmissible votes are refused by the prover */
    {
        tv_ballot b; tv_voter_secret sec;
        tv_ballot_alloc(&b, &prm); tv_secret_alloc(&sec, &prm);
        int bad[16];
        for (int a = 0; a < L; a++) bad[a] = 0;
        bad[0] = 2;
        CHECK(tv_vote(&b, &sec, NULL, &pub, keys, 5, bad, s) != 0, "prover refuses vote 2");
        if (w != TV_W_FREE) {
            for (int a = 0; a < L; a++) bad[a] = 1;
            CHECK(tv_vote(&b, &sec, NULL, &pub, keys, 5, bad, s) != 0, "prover refuses wrong weight");
        }
        tv_ballot_free(&b); tv_secret_free(&sec);
    }
    free(buf);
    free(keys);
    tv_pub_free(&pub);
}


/* a small election end to end: ballots, share checks, aggregation of every authority with
   repeated rounds until t succeed, verification, tampering, and the tally */
static void test_aggregation(int NV, c2_mode mode) {
    int n = 4, t = 3, L = 2, w = 1, d = 8;
    uint64_t q = 70368744177937ULL;
    printf("aggregation: NV=%d n=%d t=%d L=%d C2=%s\n", NV, n, t, L, mode == C2_SIGNED ? "signed" : "binary");
    tv_params prm;
    tv_params_init(&prm, n, t, L, w, d, q);
    uint8_t s[32];
    seed_from(s, 7);
    tv_pub pub;
    tv_pub_init(&pub, &prm, s);
    tv_akey keys[4];
    for (int k = 0; k < n; k++) seed_from(keys[k].key, 200 + (unsigned)k);
    size_t vn = (size_t)prm.mu * TV_N;
    poly *leaf_com = malloc(sizeof(poly) * (size_t)n * NV * prm.rows);
    poly *leaf_msg = malloc(sizeof(poly) * (size_t)n * NV * L);
    int64_t *leaf_rnd = malloc(sizeof(int64_t) * (size_t)n * NV * vn);
    long expect[8] = {0};
    int allok = 1;
    for (int i = 0; i < NV; i++) {
        tv_ballot b; tv_voter_secret sec;
        tv_ballot_alloc(&b, &prm); tv_secret_alloc(&sec, &prm);
        int v[8] = {0};
        v[(i * 7) % L] = 1;
        for (int a = 0; a < L; a++) expect[a] += v[a];
        seed_from(s, 5000 + (unsigned)i);
        tv_vote(&b, &sec, NULL, &pub, keys, (uint64_t)i, v, s);
        allok &= tv_verify_ballot(&pub, &b);
        for (int k = 1; k <= n; k++) {
            size_t idx = (size_t)(k - 1) * NV + (size_t)i;
            memcpy(leaf_com + idx * prm.rows, b.c + (size_t)k * prm.rows, sizeof(poly) * prm.rows);
            allok &= tv_check_share(&pub, &keys[k - 1], k, &b, leaf_msg + idx * L, leaf_rnd + idx * vn);
        }
        tv_ballot_free(&b); tv_secret_free(&sec);
    }
    CHECK(allok, "%d ballots valid and all shares open", NV);
    agg_params ap;
    agg_params_init(&ap, &prm, NV, mode, 64);
    agg_params_print(&ap, &prm);
    agg_contrib c[4];
    agg_state st[4];
    int succ[4] = {0}, nsucc = 0, round = 0;
    uint8_t bb1[TV_HASHBYTES];
    while (nsucc < t && round < 10) {
        round++;
        for (int k = 0; k < n; k++) {
            if (round > 1) { agg_state_free(&st[k]); agg_contrib_free(&c[k]); }
            uint8_t as[32];
            seed_from(as, 900 + (unsigned)(k + 10 * round));
            agg_stats stt;
            CHECK(agg_round1(&st[k], &c[k], &stt, &pub, &ap, k + 1, leaf_com + (size_t)k * NV * prm.rows,
                             leaf_msg + (size_t)k * NV * L, leaf_rnd + (size_t)k * NV * vn, NV, as) == 0,
                  "round %d: authority %d round 1 (%.0f ms, %zu bytes)", round, k + 1, stt.t_round1_ms, stt.bytes_round1);
        }
        agg_bb1_digest(bb1, &pub, c, n);
        nsucc = 0;
        for (int k = 0; k < n; k++) {
            agg_stats stt;
            succ[k] = agg_round2(&st[k], &c[k], &stt, &pub, &ap, bb1) == 0;
            nsucc += succ[k];
            printf("  authority %d round 2: %s after %d attempts, ||B1||/T1=%.3f ||B2||/T2=%.3f, %.0f ms, %zu bytes\n", k + 1,
                   succ[k] ? "ok" : "no attempt", stt.attempts_used, stt.shift1_ratio, stt.shift2_ratio, stt.t_round2_ms, stt.bytes_round2);
        }
    }
    CHECK(nsucc >= t, "at least t authorities succeeded after %d round(s)", round);
    int T[4], nt = 0;
    const poly *V[4];
    for (int k = 0; k < n; k++) {
        if (!succ[k]) continue;
        int ok = agg_verify(&pub, &ap, &c[k], leaf_com + (size_t)k * NV * prm.rows, NV, bb1);
        CHECK(ok, "VerAgg accepts authority %d", k + 1);
        if (nt < t) { T[nt] = k + 1; V[nt] = c[k].V; nt++; }
    }
    long counts[8];
    CHECK(tv_combine(counts, &pub, T, V, NV) == 0 && counts[0] == expect[0] && counts[L - 1] == expect[L - 1],
          "tally from authorities {%d,%d,%d}: %ld/%ld (expected %ld/%ld)", T[0], T[1], T[2], counts[0], counts[L - 1], expect[0], expect[L - 1]);
    /* tampering */
    int k0 = T[0] - 1;
    c[k0].Z1[3] += 1;
    CHECK(!agg_verify(&pub, &ap, &c[k0], leaf_com + (size_t)k0 * NV * prm.rows, NV, bb1), "reject: modified Z1");
    c[k0].Z1[3] -= 1;
    c[k0].Z2[11] -= 1;
    CHECK(!agg_verify(&pub, &ap, &c[k0], leaf_com + (size_t)k0 * NV * prm.rows, NV, bb1), "reject: modified Z2");
    c[k0].Z2[11] += 1;
    c[k0].V[0].c[0] = mod_add(c[k0].V[0].c[0], 1);
    CHECK(!agg_verify(&pub, &ap, &c[k0], leaf_com + (size_t)k0 * NV * prm.rows, NV, bb1), "reject: partial tally + 1");
    c[k0].V[0].c[0] = mod_sub(c[k0].V[0].c[0], 1);
    c[k0].inner_com[0].c[5] = mod_add(c[k0].inner_com[0].c[5], 1);
    CHECK(!agg_verify(&pub, &ap, &c[k0], leaf_com + (size_t)k0 * NV * prm.rows, NV, bb1), "reject: modified node commitment");
    c[k0].inner_com[0].c[5] = mod_sub(c[k0].inner_com[0].c[5], 1);
    c[k0].leaf_digest[0] ^= 1;
    CHECK(!agg_verify(&pub, &ap, &c[k0], leaf_com + (size_t)k0 * NV * prm.rows, NV, bb1), "reject: other set of leaves");
    c[k0].leaf_digest[0] ^= 1;
    CHECK(agg_verify(&pub, &ap, &c[k0], leaf_com + (size_t)k0 * NV * prm.rows, NV, bb1), "restored contribution verifies");
    uint8_t bb1x[TV_HASHBYTES];
    memcpy(bb1x, bb1, sizeof bb1x); bb1x[0] ^= 1;
    CHECK(!agg_verify(&pub, &ap, &c[k0], leaf_com + (size_t)k0 * NV * prm.rows, NV, bb1x), "reject: other board state");
    leaf_com[(size_t)k0 * NV * prm.rows + 3].c[0] = mod_add(leaf_com[(size_t)k0 * NV * prm.rows + 3].c[0], 1);
    CHECK(!agg_verify(&pub, &ap, &c[k0], leaf_com + (size_t)k0 * NV * prm.rows, NV, bb1), "reject: modified leaf");
    CHECK(agg_verify(&pub, &ap, &c[k0], leaf_com + (size_t)k0 * NV * prm.rows, NV, bb1) == 0, "(leaf still modified)");
    for (int k = 0; k < n; k++) { agg_state_free(&st[k]); agg_contrib_free(&c[k]); }
    free(leaf_com); free(leaf_msg); free(leaf_rnd);
    tv_pub_free(&pub);
}

static void test_evolve(void) {
    printf("EVOLVE ballots on the same core\n");
    tv_params prm;
    tv_params_init(&prm, 4, 4, 1, TV_W_FREE, 8, 70368744177937ULL);
    uint8_t s[32];
    seed_from(s, 3);
    tv_pub pub;
    tv_pub_init(&pub, &prm, s);
    ev_params ep;
    ev_params_init(&ep, &prm);
    tv_akey keys[4];
    for (int k = 0; k < 4; k++) seed_from(keys[k].key, 300 + (unsigned)k);
    size_t vn = (size_t)prm.mu * TV_N;
    poly sh[4];
    int64_t *rnd = malloc(sizeof(int64_t) * 4 * vn);
    for (int v = 0; v < 2; v++) {
        ev_ballot b;
        ev_ballot_alloc(&b, &prm);
        int att;
        seed_from(s, 40 + (unsigned)v);
        CHECK(ev_vote(&b, sh, rnd, &att, &pub, &ep, keys, 9, v, s) == 0, "EVOLVE vote %d (%d attempts)", v, att);
        CHECK(ev_verify(&pub, &ep, &b), "EVOLVE ballot %d verifies", v);
        poly acc;
        poly_zero(&acc);
        for (int j = 0; j < 4; j++) poly_add(&acc, &acc, &sh[j]);
        int ok = acc.c[0] == (uint64_t)v;
        for (int i = 1; i < TV_N; i++) ok &= acc.c[i] == 0;
        CHECK(ok, "shares add up to the vote");
        tv_sizes sz;
        size_t bytes = ev_ballot_bytes(&pub, &ep, &b, &sz);
        CHECK(bytes > 0, "EVOLVE ballot %zu bytes (com %zu, proof %zu)", bytes, sz.commitments, sz.proof);
        b.r0[1] += 1;
        CHECK(!ev_verify(&pub, &ep, &b), "reject: modified OR response");
        b.r0[1] -= 1;
        b.c[prm.d].c[0] = mod_add(b.c[prm.d].c[0], 1);
        CHECK(!ev_verify(&pub, &ep, &b), "reject: sum commits to vote + 1");
        ev_ballot_free(&b);
    }
    free(rnd);
    tv_pub_free(&pub);
}

int main(int argc, char **argv) {
    setvbuf(stdout, NULL, _IONBF, 0);
    if (argc > 1 && !strcmp(argv[1], "-q")) verbose = 0;
    test_ballots(4, 3, 1, TV_W_FREE, 8, 70368744177937ULL);
    test_ballots(5, 3, 2, 1, 8, 70368744177937ULL);
    test_ballots(3, 2, 3, 2, 7, 2199023255633ULL);
    test_evolve();
    test_aggregation(70, C2_SIGNED);
    test_aggregation(35, C2_BINARY);
    printf("%s (%d failures)\n", fails ? "FAILED" : "all protocol tests passed", fails);
    return fails != 0;
}
