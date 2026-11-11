// SJG 独立参考实现 — 核心类型与 API(设计见 ../README.md)
// C 风格 C++17:扁平数组 + CSR,无 Object/链表/层次。
#pragma once
#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

// ---------- 数据集(扁平布局) ----------
struct Dataset {
    int n = 0, m = 0;
    std::vector<int>    mtype;   // 度量类型 (0/1/4=L1族, 2=L2, 5=cos角, 3=日期, 6=编辑)
    std::vector<double> mmax;    // 归一化常数 metricmaxdis
    // 数值场: record-major 一整块 + 每记录 base + 每记录每度量 [start,end)(记录内相对下标)
    // 注意: 支持参差记录(脏数据字段缺失/错位) —— 老实现每对象自带 size 数组, 距离核用"第一个参数
    // 的下标范围同时索引双方 buffer", 此处逐位复刻该语义(含参差时的越界式读取, 落在相邻记录数据上)。
    std::vector<float>     fdata;
    std::vector<long long> rbase;  // [n] 记录 i 数值块在 fdata 的起始
    std::vector<int>       ms;     // [n*m] 度量 j 在记录 i 内的 start(相对 rbase)
    std::vector<int>       me;     // [n*m] 度量 j 在记录 i 内的 end
    // 字符串场: NUL 结尾连续池; 槽按记录内出现顺序编号
    std::vector<char>      sdata;
    std::vector<long long> soff;   // [n * nstr] 记录 i 的第 s 个字符串槽的偏移
    int nstr = 0;

    inline const char* str0(int rec) const { return sdata.data() + soff[(size_t)rec * nstr]; }   // 槽0(忠实 quirk: 编辑距离永远用槽0)
    inline const char* str_at(int rec, int slot) const { return sdata.data() + soff[(size_t)rec * nstr + slot]; }
    // cos 范数缓存: [n*m], 仅 type5 度量填充; 存"开过方并截断为float"的范数(与老核逐位一致)
    std::vector<float> cnorm;
    long long pad = 1;              // 越界零垫长度(最大记录总长; 参差记录安全网)
    std::vector<int> slen;          // [n] 槽0字符串长度缓存(编辑距离/LAZY_LB 免 strlen)
    std::vector<int> strslot;       // [m] 空间j的字符串槽号(-1=非字符串空间); 多字符串空间数据集用
};

// ---------- 查询协议(与老实现同款文件) ----------
struct Query {
    std::vector<double> w;      // 每度量权重(文件第一行)
    std::vector<int>    kvals;  // k 值列表
    std::vector<int>    pos;    // 查询 = 数据集记录下标(已 mod n)
};

// ---------- 图(CSR, 两种边模式) ----------
struct Graph {
    int n = 0, m = 0;
    std::vector<long long> off;    // [n+1]
    // 模式A(parity): per-metric 边
    std::vector<int>     enbr;
    std::vector<uint8_t> emet;
    std::vector<float>   edis;
    // 模式B(native): per-neighbor 全 m 维距离向量 (Phase C)
    bool native = false;
    std::vector<long long> noff;   // [n+1]
    std::vector<int>     nnbr;     // [U] 去重邻居
    std::vector<float>   edvec;    // [U * m] 全度量距离向量
    std::vector<std::vector<int>> level_f;   // 每度量入口表(.ind 兼容; flat ≈ {0})
    std::vector<int> perm;         // 原 id -> 内部 idx(局部性重排; 默认恒等)
    std::vector<int> rid;          // 内部 idx -> 原 id(输出用)
};

// ---------- 距离(忠实复刻 float 语义; 每调用 DC+1) ----------
struct DCCell { long long v; char pad[56]; };   // per-thread 填充计数(复刻老实现: 防伪共享)
extern DCCell g_dcp[64];
void g_dc_reset(); long long g_dc_sum();
double dist(const Dataset& D, int a, int b, int j);
bool build_graph(Dataset& D, const std::string& out_ind, const std::string& out_sjg);   // Phase B: NN-descent + skyline -> .ind (+可选 native .sjg)
bool save_native(const std::string& path, const Graph& G);   // native 二进制格式(全维距离元组为一等公民)
bool load_native(const std::string& path, const Dataset& D, Graph& G);

// ---------- 搜索 ----------
struct SearchOpt { double epsi = 1.0; bool lazy_dw = true; bool lazy_lb = false; int ef = 0; };   // ef>0 = beam模式(retset), 否则 EPSI 阈值模式
struct Result { std::vector<int> ids; std::vector<double> dists; };

struct Engine {                     // 可复用 scratch(跨查询"温热"语义与老实现一致)
    const Dataset* D = nullptr;
    const Graph*   G = nullptr;
    bool edge_bounds = false;     // EF triangle-inequality filter; explicit opt-in
    std::vector<float>   bound;     // [n*m] TI 下界, init -1
    std::vector<float>   ndist;     // [n]
    std::vector<uint8_t> vis;       // [n]
    void init(const Dataset& d, const Graph& g);
    Result search(int q, const std::vector<double>& w, int k, const SearchOpt& o);
    Result search_ef(int q, const std::vector<double>& w, int k, int ef);   // Bounded beam; optional per-edge distance bounds
};
Result brute(const Dataset& D, int q, const std::vector<double>& w, int k);

// ---------- I/O ----------
bool load_dataset(const std::string& path, Dataset& D);
bool load_query(const std::string& path, int m, Query& Q);
bool load_hjg_ind(const std::string& path, const Dataset& D, Graph& G, bool flatten_levels = false); // opt-in union of HJG levels; default remains flat-only
void locality_reorder(Dataset& D, Graph& G);   // BFS 重排节点+数据(bit-faithful: 只换内存位置)
void build_native(const Dataset& D, Graph& G);  // 模式A -> per-neighbor 全维向量(Phase C; 补算缺失度量距离)
void precompute_norms(Dataset& D);             // cos 自身范数缓存(bit-faithful)
extern bool g_use_norm, g_use_myers;
