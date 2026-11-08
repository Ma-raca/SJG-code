// 搜索引擎 — flat 联合图 best-first + TI 下界传播(逐位复刻老实现 knn_multi_soa 的 flat 路径语义)
// 语义要点(对拍关键, 详见 README):
//  - bound(=g_ansres) init -1; 种子行全列清 0; 查询末只清 w>1e-6 的列为 0(跨查询"温热")
//  - 权重阈值: 种子 1e-5 / precheck & 真距离 1e-6 / 散射 1e-4(与老实现一致, 故意不统一)
//  - float 累加语义: ndist/bound 均 float, 与 double 权重相乘后逐步截断
//  - 候选堆比较器只看 key(tie 由堆内部序决定, 与老实现一致)
//  - LAZY_DW: tau 在真距离循环前定格; LAZY_LB: 编辑长度差下界(槽0), 延后计算
#include "sjg.h"
#include <algorithm>
#include <queue>
#include <cmath>
#include <cstring>

using std::pair;
using std::vector;

void Engine::init(const Dataset& d, const Graph& g) {
    D = &d; G = &g;
    bound.assign((size_t)d.n * d.m, -1.0f);
    ndist.assign(d.n, 0.0f);
    vis.assign(d.n, 0);
}

struct CmpFirstGT {   // 只比 key —— 复刻老实现 cmpfirst_gt(tie 走堆内部序)
    bool operator()(const pair<double,int>& a, const pair<double,int>& b) const { return a.first > b.first; }
};

Result Engine::search(int qi, const vector<double>& w, int knn, const SearchOpt& o) {
    const Dataset& D = *this->D; const Graph& G = *this->G;
    const int m = D.m;
    if (knn <= 0) return {};
    std::priority_queue<pair<double,int>, vector<pair<double,int>>, CmpFirstGT> heap;
    std::priority_queue<double> distlist, que;   // max-heap: 当前 top-knn 的 distW
    // P5: 每查询预筛活跃度量(阈值与老实现各处一致: 散射1e-4 / precheck与真距离1e-6)
    int act4[256], act6[256]; int nact4 = 0, nact6 = 0; uint8_t wact4[256] = {0};
    for (int j = 0; j < m; j++) {
        if (w[j] > 0.0001)    { act4[nact4++] = j; wact4[j] = 1; }
        if (w[j] > 0.000001)  { act6[nact6++] = j; }
    }
    vector<int> tmp; tmp.reserve(4096);
    vector<int> ans; ans.reserve(256);

    // 种子(入口): 每个 w>1e-5 的度量的 level_f, 去重
    for (int j = 0; j < m; j++) if (w[j] > 0.00001)
        for (int e : G.level_f[j]) if (!vis[e]) {
            heap.push({0.0, e}); ndist[e] = 0;
            for (int k2 = 0; k2 < m; k2++) bound[(size_t)e * m + k2] = 0;
            vis[e] = 1;
        }

    while (!heap.empty()) {
        int u = heap.top().second; heap.pop();
        vis[u] = 2;   // P2: 已弹出标记 — 此后对它的下界散射是死写(行不再被读), 可跳过
        size_t n1 = (size_t)u * m;
        tmp.push_back(u);
        int flg = 0;
        // precheck: 加权下界过阈值才算真距离
        if ((int)distlist.size() < knn || ndist[u] <= distlist.top() * o.epsi) {
            float nd = 0;
            for (int t = 0; t < nact6; t++) { int j = act6[t]; nd = (float)(nd + bound[n1 + j] * w[j]); }
            ndist[u] = nd;
            if ((int)distlist.size() < knn || ndist[u] <= distlist.top() * o.epsi) flg = 1;
        }
        if (flg) {   // 真 distW(惰性早停)
            flg = 0;
            float nd = 0;
            bool pruned = false;
            double tau = (o.lazy_dw && (int)distlist.size() >= knn) ? distlist.top() * o.epsi : 1e300;
            double lb = 0; int defj[16]; int ndef = 0;
            for (int t6 = 0; t6 < nact6; t6++) {
                int j = act6[t6];
                if (o.lazy_dw && o.lazy_lb && D.mtype[j] == 6 && (D.strslot.empty() || D.strslot[j] == 0)) {   // 编辑: 长度差下界(仅槽0有slen缓存), 延后
                    int la, lb2;
                    if (!D.slen.empty()) { la = D.slen[u]; lb2 = D.slen[qi]; }
                    else { const char* s1 = D.str0(u); const char* s2 = D.str0(qi);
                           la = 0; while (s1[la]) la++; lb2 = 0; while (s2[lb2]) lb2++; }
                    lb += (double)(la > lb2 ? la - lb2 : lb2 - la) / D.mmax[j] * w[j];
                    if (ndef < 16) defj[ndef++] = j;
                    if (nd + lb > tau) { nd = (float)(nd + lb); pruned = true; break; }
                } else {
                    double dj = dist(D, u, qi, j);
                    bound[n1 + j] = (float)dj;
                    nd = (float)(nd + bound[n1 + j] * w[j]);
                    if (nd + lb > tau) { nd = (float)(nd + lb); pruned = true; break; }
                }
            }
            if (!pruned && ndef > 0) {   // 补算延后的编辑距离
                for (int di = 0; di < ndef; di++) {
                    int j = defj[di];
                    double dj = dist(D, u, qi, j);
                    bound[n1 + j] = (float)dj;
                    nd = (float)(nd + bound[n1 + j] * w[j]);
                }
                if (nd > tau) pruned = true;
            }
            ndist[u] = nd;
            if (!pruned) {
                distlist.push(nd);
                while ((int)distlist.size() > knn) distlist.pop();
                if ((int)distlist.size() < knn || nd <= distlist.top() * o.epsi) flg = 1;
            }
        }
        if (flg) {   // 扩展: TI 下界散射到邻居
            if (G.native) {   // Phase C: 每个去重邻居一次访问, 全度量同时收紧下界; 堆 key = 加权下界(更准的 frontier 排序)
                for (long long e = G.noff[u]; e < G.noff[u + 1]; e++) {
                    int nb = G.nnbr[e];
                    if (vis[nb] == 2) continue;   // P2: 已弹出 — 死写跳过
                    const float* dv = &G.edvec[(size_t)e * m];
                    size_t nn = (size_t)nb * m;
                    for (int t = 0; t < nact4; t++) {   // P5: 活跃度量列表(升序, 与逐j跳过同序)
                        int j = act4[t];
                        float diff = dv[j] - bound[n1 + j];
                        float cand = diff < 0 ? -diff : diff;
                        if (cand > bound[nn + j]) bound[nn + j] = cand;
                    }
                    if (!vis[nb]) {
                        float wb = 0;
                        for (int t = 0; t < nact4; t++) { int j = act4[t]; wb = (float)(wb + bound[nn + j] * w[j]); }
                        vis[nb] = 1; ndist[nb] = wb; heap.push({(double)wb, nb});
                    }
                }
            } else
            for (long long e = G.off[u]; e < G.off[u + 1]; e++) {
                int j = G.emet[e];
                if (!wact4[j]) continue;
                int nb = G.enbr[e];
                if (vis[nb] == 2) continue;   // P2: 已弹出 — 死写跳过
                size_t nn = (size_t)nb * m + j;
                float diff = G.edis[e] - bound[n1 + j];
                float cand = diff < 0 ? -diff : diff;
                float v = bound[nn]; if (cand > v) v = cand;
                bound[nn] = v;
                if (!vis[nb]) { vis[nb] = 1; ndist[nb] = v; heap.push({(double)v, nb}); }
            }
        }
    }
    if (distlist.empty()) { for (int u : tmp) { vis[u] = 0; for (int j = 0; j < m; j++) if (w[j] > 0.000001) bound[(size_t)u * m + j] = 0; } return {}; }   // 全零权重/空图 guard
    double minv = distlist.top();
    for (int u : tmp) {   // 收尾: 收集 ans + 复原 scratch(与老实现同款"温热"复位)
        size_t n1 = (size_t)u * m;
        vis[u] = 0;
        if (ndist[u] <= minv) { que.push(ndist[u]); ans.push_back(u); }
        for (int j = 0; j < m; j++) if (w[j] > 0.000001) bound[n1 + j] = 0;
    }
    while ((int)que.size() > knn) que.pop();
    if (que.empty()) return {};                 // NaN距离毒化比较时ans可为空; 空堆top()是UB(q_c7崩溃根因)
    double fmin = que.top();
    Result R;
    for (int u : ans) if (ndist[u] <= fmin) { R.ids.push_back(u); R.dists.push_back(ndist[u]); }
    return R;
}

Result brute(const Dataset& D, int qi, const vector<double>& w, int k) {
    std::priority_queue<pair<double,int>> top;   // max-heap, 留 k 小
    for (int i = 0; i < D.n; i++) {
        double d = 0;
        for (int j = 0; j < D.m; j++) if (w[j] > 1e-9) d += w[j] * dist(D, i, qi, j);
        if ((int)top.size() < k) top.push({d, i});
        else if (d < top.top().first) { top.pop(); top.push({d, i}); }
    }
    Result R;
    while (!top.empty()) { R.ids.push_back(top.top().second); R.dists.push_back(top.top().first); top.pop(); }
    return R;   // 降序(与老实现一致)
}

// EF(beam/retset)模式 — 经典 hnswlib 搜索语义 + native 全维下界筛(相对经典 beam 无损, 省真距离):
//  - result = 容量 ef 的最大堆(真 distW); cand = 最小堆; 终止: cand 最小 > result 最差
//  - 淘汰哲学是"名额"(算法近似, ef 是准度旋钮), 与 EPSI 的"距离证据"不同 -> 前沿由实验判
//  - 发现邻居时先用 native 全维元组算加权下界 wb: wb >= 第ef名 => 真距离也进不了名单, 直接跳过
//    (真距离 >= 下界, 与经典 beam 结果一致但省计算); 真 distW 计算带 tau=第ef名 的惰性早停
//  - visited 于发现时标记(经典语义: 被名额拒绝的节点永不复活)
Result Engine::search_ef(int qi, const vector<double>& w, int knn, int ef) {
    const Dataset& D = *this->D; const Graph& G = *this->G;
    const int m = D.m;
    if (knn <= 0 || !G.native) return {};
    if (ef < knn) ef = knn;
    int act[256]; int nact = 0;
    for (int j = 0; j < m; j++) if (w[j] > 0.000001) act[nact++] = j;
    // Keep act in metric order for seeds, edge bounds and score accumulation.
    // Only the expensive component evaluations are ordered by descending weight.
    int eval_order[256];
    for (int t = 0; t < nact; t++) eval_order[t] = t;
    std::stable_sort(eval_order, eval_order + nact,
                     [&](int a, int b) { return w[act[a]] > w[act[b]]; });
    bool reordered = false;
    for (int t = 0; t < nact; t++) reordered |= eval_order[t] != t;
    double exact_distances[256];
    std::priority_queue<pair<double,int>> result;                                            // max-heap(真distW)
    std::priority_queue<pair<double,int>, vector<pair<double,int>>, CmpFirstGT> cand;        // min-heap
    vector<int> touched; touched.reserve(4096);
    // 种子: 直接算真 distW
    for (int j = 0; j < m; j++) if (w[j] > 0.00001)
        for (int e : G.level_f[j]) if (!vis[e]) {
            vis[e] = 1; touched.push_back(e);
            size_t r = (size_t)e * m;
            double dsum = 0;
            for (int t = 0; t < nact; t++) { int j2 = act[t]; double dj = dist(D, e, qi, j2); bound[r + j2] = (float)dj; dsum += dj * w[j2]; }
            cand.push({dsum, e});
            result.push({dsum, e}); if ((int)result.size() > ef) result.pop();
        }
    while (!cand.empty()) {
        double cd = cand.top().first; int u = cand.top().second;
        if ((int)result.size() >= ef && cd > result.top().first) break;   // 经典终止
        cand.pop();
        size_t n1 = (size_t)u * m;
        for (long long e = G.noff[u]; e < G.noff[u + 1]; e++) {
            int nb = G.nnbr[e];
            if (vis[nb]) continue;
            vis[nb] = 1; touched.push_back(nb);
            const float* dv = &G.edvec[(size_t)e * m];
            double worst = ((int)result.size() >= ef) ? result.top().first : 1e300;
            double wb = 0;                                                 // 全维加权下界(免费, 来自元组+u的真距离行)
            for (int t = 0; t < nact; t++) { int j2 = act[t]; float diff = dv[j2] - bound[n1 + j2]; wb += w[j2] * (diff < 0 ? -diff : diff); }
            if (wb >= worst) continue;                                     // 下界证明进不了名单(无损筛)
            size_t r = (size_t)nb * m;
            double dsum = 0; bool ab = false;
            if (!reordered) {
                // Equal weights and already ordered weights retain the original fast path.
                for (int t = 0; t < nact; t++) {
                    int j2 = act[t];
                    double dj = dist(D, nb, qi, j2);
                    bound[r + j2] = (float)dj;
                    dsum += dj * w[j2];
                    if (dsum >= worst) { ab = true; break; }
                }
            } else {
                std::fill(exact_distances, exact_distances + nact, 0.0);
                for (int t = 0; t < nact; t++) {
                    int position = eval_order[t];
                    int j2 = act[position];
                    double dj = dist(D, nb, qi, j2);
                    exact_distances[position] = dj;
                    bound[r + j2] = (float)dj;
                    // Replay the original sum, with zero for unevaluated terms.
                    // For finite nonnegative distances and weights this is a lower
                    // bound on the original full score, without reordering its
                    // floating-point additions. Use doubles, not the float cache.
                    dsum = 0;
                    for (int s = 0; s < nact; s++)
                        dsum += exact_distances[s] * w[act[s]];
                    if (dsum >= worst) { ab = true; break; }
                }
            }
            if (ab) continue;
            // Only fully evaluated rows enter the queues and can supply edge bounds.
            cand.push({dsum, nb});
            result.push({dsum, nb}); if ((int)result.size() > ef) result.pop();
        }
    }
    for (int u : touched) vis[u] = 0;
    // 取 result 中最小的 k 个
    vector<pair<double,int>> all;
    while (!result.empty()) { all.push_back(result.top()); result.pop(); }   // 降序弹出
    Result R;
    int take = (int)all.size() < knn ? (int)all.size() : knn;
    for (int i = (int)all.size() - 1; i >= (int)all.size() - take; i--) { R.ids.push_back(all[i].second); R.dists.push_back(all[i].first); }
    return R;
}
