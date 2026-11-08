// Phase B: 并行建图 — per-metric NN-descent -> 联合 -> skyline+floor -> 直出 .ind(无链表中间态)
// 忠实复刻老实现 par_build_nndescent + skyline_prune 的确定性语义:
//  - xorshift 按 (节点,度量) 播种; 有界 top-K 池按 (dis,id) 全序, 内容与到达序无关
//  - 逆向表串行采样(升序 u, 截 rho); local join new×new + new×old, 自旋锁双向插入
//  - step3 串行发边(前向 top-emit + 反向), add_e 去重语义; level_f = [seed]+孤立点
//  - skyline: 链序候选收集(metric 升序 × 边逆插入序) -> 支配剪枝 -> Borda floor;
//    重建后最终边序 = 保留边按原插入序(老实现 prepend 两次的净效果)
// 环境变量与老实现同名: HJG_THREADS ND_K ND_EMIT ND_ITERS ND_DELTA ND_RHO HJG_SKYLINE SKY_FLOOR
// SJG_INDEX_K 仅用于索引邻居数敏感性实验; 未设置时保持生产默认值 15。
#include "sjg.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

using std::vector;

struct NDEnt { int id; double dis; bool isnew; };
static inline bool nd_less(const NDEnt& a, const NDEnt& b) {
    if (a.dis != b.dis) return a.dis < b.dis;
    return a.id < b.id;
}
static inline int nd_insert(NDEnt* L, int& cnt, int K, int cid, double d, int self) {
    if (cid == self) return 0;
    if (cnt >= K) { const NDEnt& w = L[cnt - 1]; if (d > w.dis || (d == w.dis && cid >= w.id)) return 0; }
    for (int x = 0; x < cnt; x++) if (L[x].id == cid) return 0;
    NDEnt e{cid, d, true};
    int pos = cnt; while (pos > 0 && nd_less(e, L[pos - 1])) { if (pos < K) L[pos] = L[pos - 1]; pos--; }
    if (pos < K) L[pos] = e;
    if (cnt < K) cnt++;
    return 1;
}
static inline void bl_lock(char* f) { while (__atomic_exchange_n(f, 1, __ATOMIC_ACQUIRE)) { while (__atomic_load_n(f, __ATOMIC_RELAXED)); } }
static inline void bl_unlock(char* f) { __atomic_store_n(f, 0, __ATOMIC_RELEASE); }

struct BEdge { int nbr; double dis; };

bool build_graph(Dataset& D, const std::string& out_ind, const std::string& out_sjg) {
    const int n = D.n, m = D.m;
    int max_n = 15;
    if (const char* value = getenv("SJG_INDEX_K")) {
        char* end = nullptr;
        long parsed = strtol(value, &end, 10);
        if (end == value || *end != '\0' || parsed < 1 || parsed >= n) {
            std::cerr << "[sjg build] invalid SJG_INDEX_K=" << value
                      << " (expected 1 <= K < " << n << ")\n";
            return false;
        }
        max_n = (int)parsed;
    }
    int nth  = getenv("HJG_THREADS") ? atoi(getenv("HJG_THREADS")) : m; if (nth < 1) nth = 1;
    int K    = getenv("ND_K") ? atoi(getenv("ND_K")) : (max_n + 1); if (K < max_n) K = max_n;
    int emit = getenv("ND_EMIT") ? atoi(getenv("ND_EMIT")) : max_n; if (emit > K) emit = K;
    int max_iter = getenv("ND_ITERS") ? atoi(getenv("ND_ITERS")) : 12;
    double delta = getenv("ND_DELTA") ? atof(getenv("ND_DELTA")) : 0.001;
    int rho = getenv("ND_RHO") ? atoi(getenv("ND_RHO")) : K; if (rho < 1) rho = 1; if (rho > K) rho = K;
    bool sky = getenv("HJG_SKYLINE") != NULL;
    if (getenv("SKY_ALPHA") || getenv("SKY_ALPHA_M"))
        std::cerr << "[sjg build] WARNING: SKY_ALPHA/A_M(α-RobustPrune 导航边)未在 sjg 实现, 将被忽略; 需要请用老 hjg_par\n";
    int flr  = getenv("SKY_FLOOR") ? atoi(getenv("SKY_FLOOR")) : 0;
    int cap  = getenv("SKY_CAP") ? atoi(getenv("SKY_CAP")) : 0;   // 高m时Pareto不剪(高维互不支配), Borda序截断到cap; 0=关
    std::cerr << "[sjg build] index_k=" << max_n << " nd_pool=" << K
              << " emit=" << emit << "\n";
    // SKY_MODE 消融: 先按天际线(含floor/cap)定出每节点度数预算, 再在同预算下换选边规则。
    // random = 候选中均匀随机选(排除"边少就快"); fixw = 均匀权重加权距离最小者(排除"普通固定权重剪枝就够")
    int sky_mode = 0;
    if (const char* sm = getenv("SKY_MODE")) {
        if (std::string(sm) == "random") sky_mode = 1;
        else if (std::string(sm) == "fixw") sky_mode = 2;
        else if (std::string(sm) != "skyline") { std::cerr << "[sjg build] unknown SKY_MODE=" << sm << " (skyline|random|fixw)\n"; return false; }
    }
    const int seed = 0;
    double w0 = omp_get_wtime();
    g_dc_reset();

    vector<vector<BEdge>> adj((size_t)n * m);            // adj[u*m+j] = 插入序(内存链 = 其逆)
    vector<std::unordered_set<long long>> dd(n);         // add_e 去重: key = j*n + nbr
    vector<vector<int>> lf(m);                           // level_f
    vector<vector<BEdge>> seed_emit(m);                  // 仅供入口完全无候选时兜底

    vector<NDEnt> pool((size_t)n * K);
    vector<int> cnt(n, 0);
    vector<char> flag(n, 0);
    vector<vector<int>> newl(n), oldl(n), rnew(n), rold(n);

    for (int j = 0; j < m; j++) {
        // ---- 1) 初始化: 每节点 K 个确定性随机邻居 ----
        #pragma omp parallel for schedule(dynamic, 256) num_threads(nth)
        for (int u = 0; u < n; u++) {
            NDEnt* L = &pool[(size_t)u * K]; int c = 0;
            unsigned long long s = (unsigned long long)(u + 1) * 2654435761ULL + (unsigned long long)(j + 1) * 40503ULL;
            int tries = 0;
            while (c < K && tries < K * 8 && n > 1) {
                s ^= s << 13; s ^= s >> 7; s ^= s << 17;
                int v = (int)(s % (unsigned long long)n); tries++;
                if (v == u) continue;
                bool dup = false; for (int x = 0; x < c; x++) if (L[x].id == v) { dup = true; break; }
                if (dup) continue;
                double d = dist(D, u, v, j);
                NDEnt e{v, d, true}; int pos = c; while (pos > 0 && nd_less(e, L[pos - 1])) { L[pos] = L[pos - 1]; pos--; }
                L[pos] = e; c++;
            }
            cnt[u] = c;
        }
        // ---- 2) NN-descent 迭代(增量 new/old local join) ----
        for (int it = 0; it < max_iter; it++) {
            #pragma omp parallel for schedule(dynamic, 256) num_threads(nth)
            for (int u = 0; u < n; u++) {
                newl[u].clear(); oldl[u].clear(); rnew[u].clear(); rold[u].clear();
                NDEnt* L = &pool[(size_t)u * K];
                for (int x = 0; x < cnt[u]; x++) {
                    if (L[x].isnew) { newl[u].push_back(L[x].id); L[x].isnew = false; } else oldl[u].push_back(L[x].id);
                }
            }
            for (int u = 0; u < n; u++) {                 // 逆向表: 串行升序采样(确定性)
                for (size_t x = 0; x < newl[u].size(); x++) { int v = newl[u][x]; if ((int)rnew[v].size() < rho) rnew[v].push_back(u); }
                for (size_t x = 0; x < oldl[u].size(); x++) { int v = oldl[u][x]; if ((int)rold[v].size() < rho) rold[v].push_back(u); }
            }
            long long updates = 0;
            #pragma omp parallel for schedule(dynamic, 64) num_threads(nth) reduction(+:updates)
            for (int u = 0; u < n; u++) {
                vector<int> nw = newl[u]; nw.insert(nw.end(), rnew[u].begin(), rnew[u].end());
                vector<int> od = oldl[u]; od.insert(od.end(), rold[u].begin(), rold[u].end());
                int NW = (int)nw.size(), OD = (int)od.size();
                for (int a = 0; a < NW; a++) {
                    int ia = nw[a]; NDEnt* La = &pool[(size_t)ia * K];
                    for (int b = a + 1; b < NW; b++) {
                        int ib = nw[b]; if (ia == ib) continue;
                        double d = dist(D, ia, ib, j);
                        bl_lock(&flag[ia]); int r1 = nd_insert(La, cnt[ia], K, ib, d, ia); bl_unlock(&flag[ia]);
                        bl_lock(&flag[ib]); int r2 = nd_insert(&pool[(size_t)ib * K], cnt[ib], K, ia, d, ib); bl_unlock(&flag[ib]);
                        updates += r1 + r2;
                    }
                    for (int b = 0; b < OD; b++) {
                        int ib = od[b]; if (ia == ib) continue;
                        double d = dist(D, ia, ib, j);
                        bl_lock(&flag[ia]); int r1 = nd_insert(La, cnt[ia], K, ib, d, ia); bl_unlock(&flag[ia]);
                        bl_lock(&flag[ib]); int r2 = nd_insert(&pool[(size_t)ib * K], cnt[ib], K, ia, d, ib); bl_unlock(&flag[ib]);
                        updates += r1 + r2;
                    }
                }
            }
            if (updates < (long long)(delta * (double)n * (double)K)) break;
        }
        // 保存 seed 在当前度量下已有的 top-emit 结果。正常路径不会发出这些边。
        {
            NDEnt* L = &pool[(size_t)seed * K];
            int e_ = cnt[seed] < emit ? cnt[seed] : emit;
            if (e_ < 0) e_ = 0;
            seed_emit[j].reserve(e_);
            for (int x = 0; x < e_; x++) seed_emit[j].push_back({L[x].id, L[x].dis});
        }

        // ---- 3) 发边(前向 top-emit + 反向; add_e 去重语义) + level_f ----
        lf[j].push_back(seed);
        for (int i = 0; i < n; i++) { if (i == seed) continue;
            NDEnt* L = &pool[(size_t)i * K];
            int e_ = cnt[i] < emit ? cnt[i] : emit;
            for (int x = 0; x < e_; x++) {
                int q = L[x].id; double dis = L[x].dis;
                long long kf = (long long)j * n + q, kr = (long long)j * n + i;
                if (dd[i].insert(kf).second) adj[(size_t)i * m + j].push_back({q, dis});
                if (dd[q].insert(kr).second) adj[(size_t)q * m + j].push_back({i, dis});
            }
            if (adj[(size_t)i * m + j].empty()) lf[j].push_back(i);
        }
        std::cerr << "[sjg build] metric " << j << " done  " << omp_get_wtime() - w0 << "s\n";
    }

    // 入口通常通过其他节点的反向连接获得候选。只有它在所有度量下均为空时，
    // 才发出自身已有的 top-emit 近邻，避免查询从空入口开始。
    bool seed_empty = true;
    for (int j = 0; j < m; j++)
        if (!adj[(size_t)seed * m + j].empty()) { seed_empty = false; break; }
    if (seed_empty && n > 1) {
        long long added = 0;
        for (int j = 0; j < m; j++) {
            for (const BEdge& e : seed_emit[j]) {
                int q = e.nbr;
                long long kf = (long long)j * n + q;
                long long kr = (long long)j * n + seed;
                if (dd[seed].insert(kf).second) {
                    adj[(size_t)seed * m + j].push_back(e);
                    added++;
                }
                if (dd[q].insert(kr).second) {
                    adj[(size_t)q * m + j].push_back({seed, e.dis});
                    added++;
                }
            }
        }
        if (added > 0)
            std::cerr << "[sjg build] entry fallback emitted " << added
                      << " directed candidate links\n";
    }

    // native 留存: 每节点 kept 邻居的全维距离元组(skyline 免费副产品)
    std::vector<std::vector<int>>   knbr((size_t)(out_sjg.empty() ? 0 : n));
    std::vector<std::vector<float>> kdv((size_t)(out_sjg.empty() ? 0 : n));
    const bool want_native = !out_sjg.empty();
    // ---- skyline 剪枝(逐位复刻 skyline_prune: 链序候选 + 支配 + Borda floor) ----
    if (sky) {
        double t0 = omp_get_wtime();
        long long kept_tot = 0, uni_tot = 0;
        #pragma omp parallel for schedule(dynamic) num_threads(nth) reduction(+:kept_tot,uni_tot)
        for (int u = 0; u < n; u++) {
            // 候选按链序收集: metric 升序 × 边逆插入序(= 老实现内存链走序); 哈希去重(老实现同款 unordered_map)
            vector<int> cand; vector<vector<double>> dv;
            std::unordered_map<int,int> cmap;
            for (int j = 0; j < m; j++) {
                auto& A = adj[(size_t)u * m + j];
                for (int x = (int)A.size() - 1; x >= 0; x--) {
                    int nb = A[x].nbr;
                    auto it = cmap.find(nb); int ci;
                    if (it == cmap.end()) { ci = (int)cand.size(); cmap[nb] = ci; cand.push_back(nb); dv.push_back(vector<double>(m, -1.0)); }
                    else ci = it->second;
                    dv[ci][j] = A[x].dis;
                }
            }
            int nc = (int)cand.size();
            if (nc == 0) continue;
            uni_tot += nc;
            for (int c = 0; c < nc; c++) for (int j = 0; j < m; j++)
                if (dv[c][j] < 0) dv[c][j] = dist(D, u, cand[c], j);
            vector<char> keep(nc, 0); int nkeep = 0;
            for (int a = 0; a < nc; a++) { bool dom = false;
                for (int b = 0; b < nc && !dom; b++) { if (b == a) continue; bool le = true, lt = false;
                    for (int j = 0; j < m; j++) { if (dv[b][j] > dv[a][j]) { le = false; break; }
                                                  if (dv[b][j] < dv[a][j]) lt = true; }
                    if (le && lt) dom = true; }
                if (!dom) { keep[a] = 1; nkeep++; } }
            if (flr > 0 && nkeep < flr && nkeep < nc) {       // Borda floor(与老实现同款不稳定排序语义)
                vector<int> rk(nc, 0), ord(nc);
                for (int j = 0; j < m; j++) { for (int a = 0; a < nc; a++) ord[a] = a;
                    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return dv[a][j] < dv[b][j]; });
                    for (int r = 0; r < nc; r++) rk[ord[r]] += r; }
                vector<int> rest; for (int a = 0; a < nc; a++) if (!keep[a]) rest.push_back(a);
                std::sort(rest.begin(), rest.end(), [&](int a, int b){ return rk[a] < rk[b]; });
                for (int a : rest) { if (nkeep >= flr) break; keep[a] = 1; nkeep++; } }
            if (cap > 0 && nkeep > cap) {   // SKY_CAP 截断: 保留边中 Borda 总秩最小的 cap 条
                vector<int> rk(nc, 0), ord(nc);
                for (int j = 0; j < m; j++) { for (int a = 0; a < nc; a++) ord[a] = a;
                    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return dv[a][j] < dv[b][j]; });
                    for (int r = 0; r < nc; r++) rk[ord[r]] += r; }
                vector<int> ks; for (int a = 0; a < nc; a++) if (keep[a]) ks.push_back(a);
                std::sort(ks.begin(), ks.end(), [&](int a, int b){ return rk[a] < rk[b]; });
                for (int t = cap; t < (int)ks.size(); t++) keep[ks[t]] = 0;
                nkeep = cap;
            }
            if (sky_mode != 0 && nkeep < nc) {   // 同度数预算下替换选边规则(消融)
                int budget = nkeep;
                std::fill(keep.begin(), keep.end(), 0);
                if (sky_mode == 2) {             // fixw: 均匀权重加权距离最小的 budget 个
                    vector<int> ord(nc); for (int a = 0; a < nc; a++) ord[a] = a;
                    vector<double> agg(nc, 0.0);
                    for (int a = 0; a < nc; a++) for (int j = 0; j < m; j++) {
                        double s = (j < (int)D.mmax.size() && D.mmax[j] > 0) ? D.mmax[j] : 1.0;
                        agg[a] += dv[a][j] / s;
                    }
                    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return agg[a] != agg[b] ? agg[a] < agg[b] : a < b; });
                    for (int t = 0; t < budget; t++) keep[ord[t]] = 1;
                } else {                          // random: 确定性 xorshift 按节点播种, 部分 Fisher-Yates
                    vector<int> ord(nc); for (int a = 0; a < nc; a++) ord[a] = a;
                    unsigned long long s = (unsigned long long)(u + 1) * 2654435761ULL + 97531ULL;
                    for (int t = 0; t < budget; t++) {
                        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
                        int r = t + (int)(s % (unsigned long long)(nc - t));
                        std::swap(ord[t], ord[r]); keep[ord[t]] = 1;
                    }
                }
                nkeep = budget;
            }
            kept_tot += nkeep;
            if (want_native) {   // 留存 kept 候选的全维元组(cand 序)
                for (int c = 0; c < nc; c++) if (keep[c]) {
                    knbr[u].push_back(cand[c]);
                    for (int j = 0; j < m; j++) kdv[u].push_back((float)dv[c][j]);
                }
            }
            // 重建: 最终链 = 保留边按原插入序(prepend 两次的净效果); 哈希查找
            for (int j = 0; j < m; j++) {
                auto& A = adj[(size_t)u * m + j];
                vector<BEdge> kept2; kept2.reserve(A.size());
                for (auto& e : A) {
                    auto it = cmap.find(e.nbr);
                    if (it != cmap.end() && keep[it->second]) kept2.push_back(e);
                }
                A.swap(kept2);
            }
        }
        std::cerr << "[sjg build] skyline union-deg=" << (double)uni_tot / n << " -> " << (double)kept_tot / n
                  << "  floor=" << flr << (sky_mode == 1 ? "  MODE=random" : sky_mode == 2 ? "  MODE=fixw" : "")
                  << "  " << omp_get_wtime() - t0 << "s\n";
    }

    // ---- 写 .ind(与老 writer 同构; 边序: skyline 后=原插入序, 无 skyline=逆插入序[链序]) ----
    bool io_ok = true;
    {
        std::ofstream ou(out_ind);
        if (!ou) { std::cerr << "[sjg build] cannot write " << out_ind << "\n"; return false; }
        for (int j = 0; j < m; j++) {
            ou << lf[j].size();
            for (int id : lf[j]) ou << " " << id;
            ou << "\n";
        }
        ou << n << "\n";
        for (int i = 0; i < n; i++) {
            ou << i << " " << 0 << " " << m;
            for (int j = 0; j < m; j++) {
                auto& A = adj[(size_t)i * m + j];
                if (A.empty()) { ou << " " << j << " " << 0; continue; }
                ou << " " << j << " " << 1 << " " << 0 << " " << A.size();
                if (sky) { for (auto& e : A) ou << " " << e.nbr << " " << e.dis; }
                else     { for (int x = (int)A.size() - 1; x >= 0; x--) ou << " " << A[x].nbr << " " << A[x].dis; }
            }
            ou << "\n";
        }
        ou.flush();
        if (!ou) { std::cerr << "[sjg build] write FAILED: " << out_ind << "\n"; io_ok = false; }
    }
    if (want_native) {
        double tn0 = omp_get_wtime();
        // 并行组装(产物字节不变): pass1 每节点独立生成局部缓冲(哈希去重, 去重序=读入链序),
        // pass2 前缀和 + 并行平铺拷贝。原为串行 + 线性扫(向量数据集建图慢的主因)。
        std::vector<std::vector<int>>   lnbr((size_t)n);
        std::vector<std::vector<float>> ldv((size_t)n);
        long long fills = 0;
        #pragma omp parallel for schedule(dynamic, 64) num_threads(nth) reduction(+:fills)
        for (int u = 0; u < n; u++) {
            std::unordered_map<int,int> kmap;
            for (int c = 0; c < (int)knbr[u].size(); c++) kmap[knbr[u][c]] = c;
            std::unordered_set<int> seen;
            auto& LN = lnbr[u]; auto& LD = ldv[u];
            for (int j = m - 1; j >= 0; j--) {
                auto& A = adj[(size_t)u * m + j];
                int an = (int)A.size();
                for (int t = 0; t < an; t++) {
                    int x = sky ? (an - 1 - t) : t;
                    int nb = A[x].nbr;
                    if (!seen.insert(nb).second) continue;
                    LN.push_back(nb);
                    size_t base = LD.size(); LD.resize(base + m);
                    auto it = kmap.find(nb);
                    if (it != kmap.end()) { for (int j2 = 0; j2 < m; j2++) LD[base + j2] = kdv[u][(size_t)it->second * m + j2]; }
                    else { for (int j2 = 0; j2 < m; j2++) LD[base + j2] = (float)dist(D, u, nb, j2); fills++; }   // 无skyline时兜底
                }
            }
        }
        Graph NG; NG.n = n; NG.m = m; NG.native = true;
        NG.level_f.assign(m, {});
        for (int j = 0; j < m; j++) NG.level_f[j] = lf[j];
        NG.noff.assign((size_t)n + 1, 0);
        for (int u = 0; u < n; u++) NG.noff[u + 1] = NG.noff[u] + (long long)lnbr[u].size();
        long long U = NG.noff[n];
        NG.nnbr.resize(U); NG.edvec.resize((size_t)U * m);
        #pragma omp parallel for schedule(dynamic, 256) num_threads(nth)
        for (int u = 0; u < n; u++) {
            long long b = NG.noff[u];
            for (size_t x = 0; x < lnbr[u].size(); x++) {
                NG.nnbr[b + x] = lnbr[u][x];
                for (int j2 = 0; j2 < m; j2++) NG.edvec[(size_t)(b + x) * m + j2] = ldv[u][x * m + j2];
            }
        }
        if (!save_native(out_sjg, NG)) { std::cerr << "[sjg build] native write FAILED: " << out_sjg << "\n"; io_ok = false; }
        else fprintf(stderr, "[sjg build] native -> %s (unique=%lld, fills=%lld, %.1fs)\n", out_sjg.c_str(), NG.noff[n], fills, omp_get_wtime() - tn0);
    }
    fprintf(stderr, "[sjg build] total %.1fs wall, build_dc=%lld -> %s\n",
            omp_get_wtime() - w0, g_dc_sum(), out_ind.c_str());
    return io_ok;
}
