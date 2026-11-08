// 距离核 — 逐位复刻 baselines/hjg_parallel 的 float 语义(object.cpp Object::distance)
#include "sjg.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <omp.h>

DCCell g_dcp[64] = {};
void g_dc_reset() { for (int i = 0; i < 64; i++) g_dcp[i].v = 0; }
long long g_dc_sum() { long long s = 0; for (int i = 0; i < 64; i++) s += g_dcp[i].v; return s; }
bool g_use_norm = true, g_use_myers = true;

// Myers 1999 位并行 Levenshtein(单字, 模式串<=64): 与 DP 完全同值
static int myers64(const char* P, int lp, const char* T, int lt) {
    static thread_local unsigned long long PEQ[256];
    for (int i = 0; i < lp; i++) PEQ[(unsigned char)P[i]] |= 1ULL << i;
    unsigned long long VP = ~0ULL, VN = 0, HB = 1ULL << (lp - 1);
    int score = lp;
    for (int t = 0; t < lt; t++) {
        unsigned long long EQ = PEQ[(unsigned char)T[t]];
        unsigned long long D0 = (((EQ & VP) + VP) ^ VP) | EQ | VN;
        unsigned long long HP = VN | ~(D0 | VP);
        unsigned long long HN = D0 & VP;
        if (HP & HB) score++;
        if (HN & HB) score--;
        HP = (HP << 1) | 1ULL;
        HN = HN << 1;
        VP = HN | ~(D0 | HP);
        VN = D0 & HP;
    }
    for (int i = 0; i < lp; i++) PEQ[(unsigned char)P[i]] = 0;   // 只清用过的项(避免 2KB memset)
    return score;
}

// 滚动 DP(老实现同款), 供 >64 回退与自测
int edit_dp(const char* s1, int len1, const char* s2, int len2) {
    static thread_local std::vector<int> prev, cur;
    if ((int)prev.size() < len2 + 1) { prev.resize(len2 + 1); cur.resize(len2 + 1); }
    for (int jj = 0; jj <= len2; jj++) prev[jj] = jj;
    for (int i = 1; i <= len1; i++) {
        cur[0] = i;
        const char ci = s1[i - 1];
        for (int jj = 1; jj <= len2; jj++) {
            int cost = (ci == s2[jj - 1]) ? 0 : 1;
            int v = prev[jj] + 1;
            int v2 = cur[jj - 1] + 1; if (v2 < v) v = v2;
            int v3 = prev[jj - 1] + cost; if (v3 < v) v = v3;
            cur[jj] = v;
        }
        std::swap(prev, cur);
    }
    return prev[len2];
}
int edit_myers(const char* s1, int len1, const char* s2, int len2) {   // 精确同值(Levenshtein 对称)
    if (len1 == 0) return len2;
    if (len2 == 0) return len1;
    if (len1 <= 64) return myers64(s1, len1, s2, len2);
    if (len2 <= 64) return myers64(s2, len2, s1, len1);
    return edit_dp(s1, len1, s2, len2);
}

static double trans_(double ddd1, double ddd2) {
    int dd1 = (int)ddd1, dd2 = (int)ddd2;
    int m1 = ((dd1 % 10000) / 100 + 9) % 12;
    int y1 = dd1 / 10000 - m1 / 10;
    int d1 = 365 * y1 + y1 / 4 - y1 / 100 + y1 / 400 + (m1 * 306 + 5) / 10 + (dd1 % 100 - 1);
    int m2 = ((dd2 % 10000) / 100 + 9) % 12;
    int y2 = (dd2 / 10000) - m2 / 10;
    int d2 = 365 * y2 + y2 / 4 - y2 / 100 + y2 / 400 + (m2 * 306 + 5) / 10 + (dd2 % 100 - 1);
    return d1 - d2;
}
static double trans2_(double ddd1, double ddd2) {
    int dd1 = (int)ddd1, dd2 = (int)ddd2;
    int y1 = dd1 / 10000, m1 = (dd1 / 100) % 100, d1 = dd1 % 100;
    int y2 = dd2 / 10000, m2 = (dd2 / 100) % 100, d2 = dd2 % 100;
    return ((y1 - y2) * 60 + m1 - m2) * 60 + d1 - d2;
}

static double jaccard_cstr(const char* a, int la, const char* b, int lb) {
    // FNV-1a 哈希令牌 + 小数组排序归并; 集合小(1~20令牌), 每次调用现算比 edit 还便宜
    auto toks = [](const char* s, int len, uint64_t* out, int cap) {
        int n = 0; uint64_t h = 1469598103934665603ULL; bool in_tok = false;
        for (int i = 0; i < len; i++) {
            char c = s[i];
            if (c == '|') { if (in_tok && n < cap) out[n++] = h; h = 1469598103934665603ULL; in_tok = false; }
            else { h = (h ^ (uint64_t)(unsigned char)c) * 1099511628211ULL; in_tok = true; }
        }
        if (in_tok && n < cap) out[n++] = h;
        std::sort(out, out + n);
        n = (int)(std::unique(out, out + n) - out);
        return n;
    };
    constexpr int CAP = 256;
    uint64_t ha[CAP], hb[CAP];
    int na = toks(a, la, ha, CAP), nb = toks(b, lb, hb, CAP);
    if (na == 0 && nb == 0) return 0.0;
    int i = 0, j2 = 0, inter = 0;
    while (i < na && j2 < nb) {
        if (ha[i] == hb[j2]) { inter++; i++; j2++; }
        else if (ha[i] < hb[j2]) i++;
        else j2++;
    }
    return 1.0 - (double)inter / (double)(na + nb - inter);
}

double dist(const Dataset& D, int a, int b, int j) {
    g_dcp[omp_get_thread_num() & 63].v++;
    double dis = 0;
    switch (D.mtype[j]) {
    case 0: case 1: case 4: {                       // L1 族: float 累加, 无分支 abs
        const float* pa = D.fdata.data() + D.rbase[a];   // 老语义: 用 a 的下标范围索引双方 buffer
        const float* pb = D.fdata.data() + D.rbase[b];
        int s0 = D.ms[(size_t)a * D.m + j], e0 = D.me[(size_t)a * D.m + j];
        float fdis = 0;
        for (int i = s0; i < e0; i++) { float diff = pa[i] - pb[i]; fdis += diff < 0 ? -diff : diff; }
        dis = fdis;
        break; }
    case 5: {                                       // cos -> 角度(度): float 累加
        const float* pa = D.fdata.data() + D.rbase[a];
        const float* pb = D.fdata.data() + D.rbase[b];
        int s0 = D.ms[(size_t)a * D.m + j], e0 = D.me[(size_t)a * D.m + j];
        // B2 忠实性: 老核把 pow 结果截断回 float、除法在 float 精度做(sa1/sa2/sa3 全 float)
        float n1, n2, sa3 = 0;
        if (g_use_norm && !D.cnorm.empty()
            && D.ms[(size_t)b * D.m + j] == s0 && D.me[(size_t)b * D.m + j] == e0) {   // 范数缓存(已开方截断): 双方范围一致才可用
            n1 = D.cnorm[(size_t)a * D.m + j]; n2 = D.cnorm[(size_t)b * D.m + j];
            for (int i = s0; i < e0; i++) { sa3 += pa[i] * pb[i]; }
        } else {
            float sa1 = 0, sa2 = 0;
            for (int i = s0; i < e0; i++) { float x = pa[i], y = pb[i]; sa1 += x * x; sa2 += y * y; sa3 += x * y; }
            n1 = (float)pow((double)sa1, 0.5); n2 = (float)pow((double)sa2, 0.5);
        }
        float prod = n1 * n2;
        if (prod == 0) { dis = (n1 == n2) ? 0.0 : 180.0; break; }   // 零向量守卫: 同零=0, 零对非零=最大角
                                                 // (老实现打印Error!!!后照除, 0/0=NaN毒化比较; 有意偏离, 见文档)
        dis = sa3 / prod;                       // float 除法(老语义), 再宽化为 double
        if (dis > 1) dis = 0.99999999999999999;
        dis = fabs(acos(dis) * 180 / 3.1415926);
        break; }
    case 2: {                                       // L2: float 平方累加, pow(.,0.5)
        const float* pa = D.fdata.data() + D.rbase[a];
        const float* pb = D.fdata.data() + D.rbase[b];
        int s0 = D.ms[(size_t)a * D.m + j], e0 = D.me[(size_t)a * D.m + j];
        float fdis = 0;
        for (int i = s0; i < e0; i++) { float diff = pa[i] - pb[i]; fdis += diff * diff; }
        dis = fdis; dis = pow(dis, 0.5);
        break; }
    case 3: {                                       // 日期(+时间)
        const float* pa = D.fdata.data() + D.rbase[a];
        const float* pb = D.fdata.data() + D.rbase[b];
        int s0 = D.ms[(size_t)a * D.m + j], e0 = D.me[(size_t)a * D.m + j];
        dis += trans_(pa[s0], pb[s0]);
        if (e0 - s0 > 1) dis = dis * 24 * 3600 + trans2_(pa[s0 + 1], pb[s0 + 1]);
        dis = fabs(dis);
        break; }
    case 6: {                                       // 编辑距离(单字符串空间=槽0, 与忠实quirk一致; 多字符串空间取本空间槽位)
        int slot = D.strslot.empty() ? 0 : D.strslot[j];
        const char* s1 = D.str_at(a, slot); const char* s2 = D.str_at(b, slot);
        int len1, len2;
        if (slot == 0 && !D.slen.empty()) { len1 = D.slen[a]; len2 = D.slen[b]; }
        else { len1 = 0; while (s1[len1]) len1++; len2 = 0; while (s2[len2]) len2++; }
        if (len1 == 0) return len2;
        if (len2 == 0) return len1;
        dis = g_use_myers ? edit_myers(s1, len1, s2, len2) : edit_dp(s1, len1, s2, len2);
        break; }
    case 7: {                                       // Jaccard 令牌集(串存 "a|b|c"; 1-|∩|/|∪|; 取本空间槽位)
        int slot = D.strslot.empty() ? 0 : D.strslot[j];
        const char* s1 = D.str_at(a, slot); const char* s2 = D.str_at(b, slot);
        int len1, len2;
        if (slot == 0 && !D.slen.empty()) { len1 = D.slen[a]; len2 = D.slen[b]; }
        else { len1 = 0; while (s1[len1]) len1++; len2 = 0; while (s2[len2]) len2++; }
        dis = jaccard_cstr(s1, len1, s2, len2);
        break; }
    default:
        std::cout << "Metic Fail " << std::endl;
        break;
    }
    return dis * 1.0 / D.mmax[j];
}
