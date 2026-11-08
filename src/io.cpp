// I/O — 数据集(data_2med.txt 同款)、查询协议、hjg .ind 导入(复刻读入链序)
#include "sjg.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <cstdlib>

using std::string;
using std::vector;

bool load_dataset(const string& path, Dataset& D) {
    std::ifstream in(path);
    if (!in) { std::cerr << "[sjg] cannot open data: " << path << "\n"; return false; }
    in >> D.n >> D.m;
    if (!in || D.n <= 0 || D.m <= 0 || D.m > 255) { std::cerr << "[sjg] bad header n=" << D.n << " m=" << D.m << "\n"; return false; }
    D.mtype.resize(D.m); D.mmax.resize(D.m);
    for (int j = 0; j < D.m; j++) in >> D.mtype[j];
    for (int j = 0; j < D.m; j++) in >> D.mmax[j];
    string line; getline(in, line);

    vector<int> strslot(D.m, -1);
    D.nstr = 0;
    for (int j = 0; j < D.m; j++) if (D.mtype[j] == 6 || D.mtype[j] == 7) strslot[j] = D.nstr++;
    D.strslot = strslot;
    D.soff.assign((size_t)D.n * std::max(D.nstr, 1), 0);
    if (D.nstr > 0) D.slen.assign(D.n, 0);
    D.rbase.resize(D.n);
    D.ms.resize((size_t)D.n * D.m);
    D.me.resize((size_t)D.n * D.m);

    long long maxrec = 1;   // B1: 垫长 = 最大记录总长(参差记录用第一参数范围索引双方, 需覆盖整条记录)
    for (int i = 0; i < D.n; i++) {
        if (!getline(in, line)) { std::cerr << "[sjg] premature EOF at record " << i << "\n"; return false; }
        // P4: strtod 游标解析(imagenet 42.6s -> 秒级); 语义严格复刻 istream:
        //  - 缺失尾字段 = 空 token(字符串=空串, 数值=单个0)
        //  - 数值 token 逗号当空白, 逐值解析, 一旦失败粘性失败(其后全 0, 复刻 failbit)
        //  - istream 拒绝 inf/nan 文本 -> 同样按失败处理
        const char* p = line.c_str(); const char* lend = p + line.size();
        static thread_local string tok;
        D.rbase[i] = (long long)D.fdata.size();
        int fpos = 0;
        for (int j = 0; j < D.m; j++) {
            while (p < lend && (*p == ' ' || *p == '\t' || *p == '\r')) p++;
            const char* ts = p;
            while (p < lend && *p != ' ' && *p != '\t' && *p != '\r') p++;
            tok.assign(ts, p);   // 空 token(缺失字段)= ""
            if (D.mtype[j] == 6 || D.mtype[j] == 7) {
                D.soff[(size_t)i * D.nstr + strslot[j]] = (long long)D.sdata.size();
                for (char c : tok) D.sdata.push_back(c);
                D.sdata.push_back('\0');
                if (strslot[j] == 0) D.slen[i] = (int)tok.size();
                D.ms[(size_t)i * D.m + j] = fpos; D.me[(size_t)i * D.m + j] = fpos;
            } else {
                int dim = 1;
                for (char& c : tok) if (c == ',') { c = ' '; dim++; }
                D.ms[(size_t)i * D.m + j] = fpos;
                const char* q = tok.c_str();
                bool failed = false;
                for (int k = 0; k < dim; k++) {
                    double v = 0;
                    if (!failed) {
                        const char* r = q; while (*r == ' ') r++;
                        if (*r == 'n' || *r == 'N' || *r == 'i' || *r == 'I') failed = true;   // istream 拒 inf/nan
                        else { char* e2; v = strtod(q, &e2); if (e2 == q || e2 == r) { v = 0; failed = true; } else q = e2; }
                        if (failed) v = 0;
                    }
                    D.fdata.push_back((float)v);
                }
                fpos += dim;
                D.me[(size_t)i * D.m + j] = fpos;
            }
        }
        if (fpos > maxrec) maxrec = fpos;
    }
    // 参差记录(脏数据)会引发老实现式的越界读: 末尾补零垫(最大记录总长), 保证最后一条记录也读不出数组
    D.pad = maxrec;
    for (long long k = 0; k < maxrec; k++) D.fdata.push_back(0.0f);
    return true;
}

bool load_query(const string& path, int m, Query& Q) {
    std::ifstream in(path);
    if (!in) { std::cerr << "[sjg] cannot open query: " << path << "\n"; return false; }
    Q.w.resize(m);
    for (int j = 0; j < m; j++) in >> Q.w[j];
    int quemk; in >> quemk;
    Q.kvals.resize(quemk);
    for (int i = 0; i < quemk; i++) in >> Q.kvals[i];
    int quemr; in >> quemr;
    double notd;
    for (int i = 0; i < quemr; i++) in >> notd;
    int qcount; in >> qcount;
    if (!in || qcount <= 0 || qcount > 1000000) { std::cerr << "[sjg] bad query qcount\n"; return false; }
    Q.pos.resize(qcount);
    for (int i = 0; i < qcount; i++) in >> Q.pos[i];
    if (!in) { std::cerr << "[sjg] truncated query file\n"; return false; }
    for (int k : Q.kvals) if (k <= 0) { std::cerr << "[sjg] bad k\n"; return false; }
    for (int& p : Q.pos) { p %= (1 << 30); if (p < 0) { std::cerr << "[sjg] negative query pos\n"; return false; } }
    return true;
}

// hjg .ind(flat)导入 -> 模式A CSR。
// 老实现读 .ind 用前插链表(metric 前插、level 前插、edge 前插), build_soa 再按链序拍平;
// 等价于: metric 逆文件序 × (level 逆文件序) × edge 逆文件序 —— 此处逐位复刻该边序。
bool load_hjg_ind(const string& path, const Dataset& D, Graph& G, bool flatten_levels) {
    std::ifstream in(path);
    if (!in) { std::cerr << "[sjg] cannot open ind: " << path << "\n"; return false; }
    G.n = D.n; G.m = D.m; G.native = false;
    G.level_f.assign(D.m, {});
    for (int j = 0; j < D.m; j++) {
        int num; in >> num;
        if (!in || num < 0 || num > D.n) return false;
        for (int t = 0; t < num; t++) {
            int id; in >> id;
            if (!in || id < 0 || id >= D.n) return false;
            G.level_f[j].push_back(id);
        }
    }
    int n2; in >> n2;
    if (n2 != D.n) { std::cerr << "[sjg] ind n=" << n2 << " != data n=" << D.n << "\n"; return false; }
    G.off.assign((size_t)D.n + 1, 0);
    // 两遍太麻烦: 先收进临时 per-node 缓冲(小), 再拼 CSR
    struct E { int nbr; uint8_t met; float dis; };
    vector<vector<E>> tmp((size_t)D.n);
    struct MetBlock { int met; vector<vector<E>> levels; };
    vector<char> seen(D.n, 0);
    long long nonzero_node_levels = 0, nonzero_edge_blocks = 0;
    for (int i = 0; i < D.n; i++) {
        int id, level, ne;
        in >> id >> level >> ne;
        if (!in || id < 0 || id >= D.n || ne < 0 || ne > D.m || level < 0 || seen[id]) { std::cerr << "[sjg] corrupt .ind at node block " << i << "\n"; return false; }
        seen[id] = 1;
        if (level != 0) {
            if (!flatten_levels) { std::cerr << "[sjg] non-flat .ind (node " << id << " level " << level << ") — refuse\n"; return false; }
            nonzero_node_levels++;
        }
        vector<MetBlock> blocks(ne);
        for (int b = 0; b < ne; b++) {
            int nl; in >> blocks[b].met >> nl;
            if (!in || blocks[b].met < 0 || blocks[b].met >= D.m || nl < 0 || nl > 4096) { std::cerr << "[sjg] corrupt .ind metric block\n"; return false; }
            blocks[b].levels.resize(nl);
            for (int l = 0; l < nl; l++) {
                int lv, nd; in >> lv >> nd;
                if (!in || lv < 0 || nd < 0 || nd > D.n) { std::cerr << "[sjg] corrupt .ind level block\n"; return false; }
                if (lv != 0) {
                    if (!flatten_levels) { std::cerr << "[sjg] non-flat edge level — refuse\n"; return false; }
                    nonzero_edge_blocks++;
                }
                auto& es = blocks[b].levels[l];
                es.resize(nd);
                for (int e = 0; e < nd; e++) {
                    int id2; std::string dt; in >> id2 >> dt;
                    if (!in || id2 < 0 || id2 >= D.n) { std::cerr << "[sjg] corrupt .ind edge (node " << id << ")\n"; return false; }
                    double dd = strtod(dt.c_str(), nullptr);   // 接受 "nan"/"inf"(脏记录的 NaN 距离; istream>> 会拒绝)
                    es[e] = { id2, (uint8_t)blocks[b].met, (float)dd };
                }
            }
        }
        // 复刻链序: metric 逆序 × level 逆序 × edge 逆序
        auto& out = tmp[id];
        for (int b = ne - 1; b >= 0; b--)
            for (int l = (int)blocks[b].levels.size() - 1; l >= 0; l--) {
                auto& es = blocks[b].levels[l];
                for (int e = (int)es.size() - 1; e >= 0; e--) { es[e].met = (uint8_t)blocks[b].met; out.push_back(es[e]); }
            }
    }
    long long tot = 0;
    for (int i = 0; i < D.n; i++) { G.off[i] = tot; tot += (long long)tmp[i].size(); }
    G.off[D.n] = tot;
    G.enbr.resize(tot); G.emet.resize(tot); G.edis.resize(tot);
    for (int i = 0; i < D.n; i++) {
        long long base = G.off[i];
        for (size_t e = 0; e < tmp[i].size(); e++) {
            G.enbr[base + e] = tmp[i][e].nbr;
            G.emet[base + e] = tmp[i][e].met;
            G.edis[base + e] = tmp[i][e].dis;
        }
    }
    G.perm.resize(D.n); G.rid.resize(D.n);
    for (int i = 0; i < D.n; i++) { G.perm[i] = i; G.rid[i] = i; }
    std::cerr << "[sjg] ind loaded: n=" << D.n << " edges=" << tot << " avg_deg=" << (double)tot / D.n << "\n";
    if (flatten_levels) std::cerr << "[sjg] HJG_FLATTEN_LEVELS: union of all stored levels; nonzero node levels="
        << nonzero_node_levels << " nonzero edge blocks=" << nonzero_edge_blocks
        << "; original entry lists retained; no Skyline pruning\n";
    return true;
}

// BFS 局部性重排: 只重排"内存位置"(节点状态/边/数据), 不改任何数值与遍历序 -> 结果逐位不变。
// 邻居在内存中相邻 -> bound/ndist/vis 散射与数据读取的 cache/TLB 局部性变好(老实现 build_soa 同思路)。
void locality_reorder(Dataset& D, Graph& G) {
    const int n = D.n, m = D.m;
    std::vector<int> order; order.reserve(n);
    std::vector<char> seen((size_t)n, 0);
    int start = 0;
    for (auto& lf : G.level_f) if (!lf.empty()) { start = lf[0]; break; }
    order.push_back(start); seen[start] = 1;
    for (size_t h = 0; h < order.size(); h++) {
        int u = order[h];
        if (G.native) {
            for (long long e = G.noff[u]; e < G.noff[u + 1]; e++) {
                int nb = G.nnbr[e];
                if (!seen[nb]) { seen[nb] = 1; order.push_back(nb); }
            }
        } else
        for (long long e = G.off[u]; e < G.off[u + 1]; e++) {
            int nb = G.enbr[e];
            if (!seen[nb]) { seen[nb] = 1; order.push_back(nb); }
        }
    }
    for (int i = 0; i < n; i++) if (!seen[i]) order.push_back(i);
    std::vector<int> perm((size_t)n);          // 原 id -> 新 idx
    for (int idx = 0; idx < n; idx++) perm[order[idx]] = idx;
    // --- 重排图(边序逐节点保持不变; 邻居 id 映射到新 idx) ---
    if (G.native) {
        std::vector<long long> noff2((size_t)n + 1, 0);
        std::vector<int> nnbr2(G.nnbr.size());
        std::vector<float> edvec2(G.edvec.size());
        long long p = 0;
        for (int idx = 0; idx < n; idx++) {
            int old = order[idx]; noff2[idx] = p;
            for (long long e = G.noff[old]; e < G.noff[old + 1]; e++) {
                nnbr2[p] = perm[G.nnbr[e]];
                for (int j = 0; j < m; j++) edvec2[(size_t)p * m + j] = G.edvec[(size_t)e * m + j];
                p++;
            }
        }
        noff2[n] = p;
        G.noff.swap(noff2); G.nnbr.swap(nnbr2); G.edvec.swap(edvec2);
    } else {
        std::vector<long long> off2((size_t)n + 1, 0);
        std::vector<int> enbr2(G.enbr.size());
        std::vector<uint8_t> emet2(G.emet.size());
        std::vector<float> edis2(G.edis.size());
        long long p = 0;
        for (int idx = 0; idx < n; idx++) {
            int old = order[idx]; off2[idx] = p;
            for (long long e = G.off[old]; e < G.off[old + 1]; e++) {
                enbr2[p] = perm[G.enbr[e]]; emet2[p] = G.emet[e]; edis2[p] = G.edis[e]; p++;
            }
        }
        off2[n] = p;
        G.off.swap(off2); G.enbr.swap(enbr2); G.emet.swap(emet2); G.edis.swap(edis2);
    }
    for (auto& lf : G.level_f) for (int& x : lf) x = perm[x];
    // --- 物理重排数据(数值块 + 字符串池): 值不变, 只换位置 ---
    {
        std::vector<float> fdata2; fdata2.reserve(D.fdata.size());
        std::vector<long long> rbase2((size_t)n);
        std::vector<int> ms2((size_t)n * m), me2((size_t)n * m);
        for (int idx = 0; idx < n; idx++) {
            int old = order[idx];
            rbase2[idx] = (long long)fdata2.size();
            long long b = D.rbase[old];
            long long len = 0;
            for (int j = 0; j < m; j++) { ms2[(size_t)idx * m + j] = D.ms[(size_t)old * m + j];
                                          me2[(size_t)idx * m + j] = D.me[(size_t)old * m + j];
                                          if (D.me[(size_t)old * m + j] > len) len = D.me[(size_t)old * m + j]; }
            for (long long t = 0; t < len; t++) fdata2.push_back(D.fdata[b + t]);
        }
        for (long long t = 0; t < D.pad; t++) fdata2.push_back(0.0f);   // 越界垫: 与 load 同源(最大记录总长)
        D.fdata.swap(fdata2); D.rbase.swap(rbase2); D.ms.swap(ms2); D.me.swap(me2);
        if (D.nstr > 0) {
            std::vector<char> sdata2; sdata2.reserve(D.sdata.size());
            std::vector<long long> soff2((size_t)n * D.nstr);
            std::vector<int> slen2((size_t)n, 0);
            for (int idx = 0; idx < n; idx++) {
                int old = order[idx];
                for (int s = 0; s < D.nstr; s++) {
                    soff2[(size_t)idx * D.nstr + s] = (long long)sdata2.size();
                    const char* p = D.sdata.data() + D.soff[(size_t)old * D.nstr + s];
                    int L = 0;
                    while (*p) { sdata2.push_back(*p++); L++; }
                    sdata2.push_back('\0');
                    if (s == 0) slen2[idx] = L;
                }
            }
            D.sdata.swap(sdata2); D.soff.swap(soff2); D.slen.swap(slen2);
        }
    }
    G.perm.swap(perm);
    G.rid = order;
    std::cerr << "[sjg] locality reorder done (BFS from " << start << ")\n";
}

// cos(type5)自身范数缓存: 与核内同序 float 累加 -> 值逐位一致
void precompute_norms(Dataset& D) {
    bool any = false;
    for (int j = 0; j < D.m; j++) if (D.mtype[j] == 5) any = true;
    if (!any) return;
    D.cnorm.assign((size_t)D.n * D.m, 0.0f);
    for (int i = 0; i < D.n; i++) {
        const float* pa = D.fdata.data() + D.rbase[i];
        for (int j = 0; j < D.m; j++) if (D.mtype[j] == 5) {
            int s0 = D.ms[(size_t)i * D.m + j], e0 = D.me[(size_t)i * D.m + j];
            float sa = 0;
            for (int t = s0; t < e0; t++) { float x = pa[t]; sa += x * x; }
            D.cnorm[(size_t)i * D.m + j] = (float)pow((double)sa, 0.5);   // 开方并截断(与老核 sa1=pow(sa1,0.5) 同语义)
        }
    }
    std::cerr << "[sjg] cos norms cached\n";
}

// Phase C: 模式A(per-metric 边) -> native(per-neighbor 全 m 维距离向量)。
// skyline 建图本已算过全向量但被 .ind 丢弃; 此处补算缺失度量(一次性索引成本)。
// 收益机制: TI 下界在所有度量上同时传播 -> 下界更紧 -> precheck 拦截更多真距离计算(打 DC)。
void build_native(const Dataset& D, Graph& G) {
    const int n = D.n, m = D.m;
    G.noff.assign((size_t)n + 1, 0);
    std::vector<int> nnbr; std::vector<float> dvec;
    nnbr.reserve(G.enbr.size()); dvec.reserve(G.enbr.size() * 2);
    long long computed = 0;
    std::unordered_map<int, size_t> uniq;   // P1: 单遍哈希去重(原 O(E×U) 双重扫描, hub 图被压垮)
    for (int u = 0; u < n; u++) {
        G.noff[u] = (long long)nnbr.size();
        uniq.clear();
        for (long long t = G.off[u], e = G.off[u + 1]; t < e; t++) {
            int nb = G.enbr[t];
            size_t x;
            auto it = uniq.find(nb);
            if (it == uniq.end()) {
                x = nnbr.size(); uniq[nb] = x; nnbr.push_back(nb);
                dvec.resize(dvec.size() + m, -1.0f);
            } else x = it->second;
            dvec[x * m + G.emet[t]] = G.edis[t];   // add_e 去重保证每(metric,nbr)至多一条 -> 无冲突
        }
        for (long long x = G.noff[u]; x < (long long)nnbr.size(); x++)
            for (int j = 0; j < m; j++) if (dvec[(size_t)x * m + j] < 0) {
                dvec[(size_t)x * m + j] = (float)dist(D, u, nnbr[x], j); computed++;
            }
    }
    G.noff[n] = (long long)nnbr.size();
    G.nnbr.swap(nnbr); G.edvec.swap(dvec);
    G.native = true;
    std::cerr << "[sjg] native: unique_nbrs=" << G.nnbr.size() << " (per-metric edges=" << G.enbr.size()
              << ", dedup " << (double)G.enbr.size() / (double)G.nnbr.size() << "x), filled=" << computed << "\n";
}

// ---- native 二进制格式 .sjg ----
// magic "SJG1" | int n, m | 每度量: int cnt + ids (level_f) | int64 noff[n+1] | int32 nnbr[U] | float edvec[U*m]
// 二进制浮点: NaN(脏记录)天然无损, 无 istream 文本坑。
bool save_native(const std::string& path, const Graph& G) {
    std::ofstream ou(path, std::ios::binary);
    if (!ou) { std::cerr << "[sjg] cannot write " << path << "\n"; return false; }
    ou.write("SJG1", 4);
    ou.write((const char*)&G.n, 4); ou.write((const char*)&G.m, 4);
    for (int j = 0; j < G.m; j++) {
        int c = (int)G.level_f[j].size();
        ou.write((const char*)&c, 4);
        if (c) ou.write((const char*)G.level_f[j].data(), 4LL * c);
    }
    ou.write((const char*)G.noff.data(), 8LL * (G.n + 1));
    long long U = G.noff[G.n];
    ou.write((const char*)G.nnbr.data(), 4LL * U);
    ou.write((const char*)G.edvec.data(), 4LL * U * G.m);
    return (bool)ou;
}
bool load_native(const std::string& path, const Dataset& D, Graph& G) {
    std::ifstream in(path, std::ios::binary);
    if (!in) { std::cerr << "[sjg] cannot open " << path << "\n"; return false; }
    char magic[4] = {0}; in.read(magic, 4);
    if (!in || std::string(magic, 4) != "SJG1") { std::cerr << "[sjg] bad magic\n"; return false; }
    in.read((char*)&G.n, 4); in.read((char*)&G.m, 4);
    if (!in || G.n != D.n || G.m != D.m) { std::cerr << "[sjg] .sjg n/m mismatch\n"; return false; }
    G.level_f.assign(G.m, {});
    for (int j = 0; j < G.m; j++) {
        int c = -1; in.read((char*)&c, 4);
        if (!in || c < 0 || c > G.n) { std::cerr << "[sjg] corrupt .sjg level_f\n"; return false; }
        G.level_f[j].resize(c);
        if (c) in.read((char*)G.level_f[j].data(), 4LL * c);
        for (int x : G.level_f[j]) if (x < 0 || x >= G.n) { std::cerr << "[sjg] corrupt .sjg entry id\n"; return false; }
    }
    G.noff.resize((size_t)G.n + 1);
    in.read((char*)G.noff.data(), 8LL * (G.n + 1));
    if (!in) { std::cerr << "[sjg] truncated .sjg (noff)\n"; return false; }
    long long U = G.noff[G.n];
    if (U < 0 || U > (long long)4e9) { std::cerr << "[sjg] corrupt .sjg U=" << U << "\n"; return false; }
    for (int i = 0; i < G.n; i++) if (G.noff[i] < 0 || G.noff[i] > G.noff[i + 1]) { std::cerr << "[sjg] corrupt .sjg noff\n"; return false; }
    G.nnbr.resize(U); G.edvec.resize(U * (size_t)G.m);
    in.read((char*)G.nnbr.data(), 4LL * U);
    in.read((char*)G.edvec.data(), 4LL * U * G.m);
    if (!in) { std::cerr << "[sjg] truncated .sjg (edges)\n"; return false; }
    for (long long t = 0; t < U; t++) if (G.nnbr[t] < 0 || G.nnbr[t] >= G.n) { std::cerr << "[sjg] corrupt .sjg nbr id\n"; return false; }
    G.native = true;
    G.perm.resize(G.n); G.rid.resize(G.n);
    for (int i = 0; i < G.n; i++) { G.perm[i] = i; G.rid[i] = i; }
    std::cerr << "[sjg] native loaded: n=" << G.n << " unique_nbrs=" << U << " avg=" << (double)U / G.n << "\n";
    return (bool)in;
}
