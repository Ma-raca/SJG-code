// CLI + 查询循环 — 输出协议与老实现 ans-tbg.txt 逐行同款(k/D/t/dc), 现有评测脚本零改动
#include "sjg.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <omp.h>


static int external_nq=0;
static void external_append(Dataset& D) {
 const char* path=getenv("EXTERNAL_QUERY_DATA"); if(!path)return;
 Dataset X; if(!load_dataset(path,X)||D.m!=X.m||D.mtype!=X.mtype||D.mmax!=X.mmax) {fprintf(stderr,"external metadata mismatch\n");exit(2);}
 precompute_norms(X); external_nq=X.n;
 long long fb=D.fdata.size(),sb=D.sdata.size();
 D.fdata.insert(D.fdata.end(),X.fdata.begin(),X.fdata.end());
 D.sdata.insert(D.sdata.end(),X.sdata.begin(),X.sdata.end());
 for(auto v:X.rbase)D.rbase.push_back(fb+v);
 for(auto v:X.soff)D.soff.push_back(sb+v);
 D.ms.insert(D.ms.end(),X.ms.begin(),X.ms.end());D.me.insert(D.me.end(),X.me.begin(),X.me.end());
 D.slen.insert(D.slen.end(),X.slen.begin(),X.slen.end());D.cnorm.insert(D.cnorm.end(),X.cnorm.begin(),X.cnorm.end());
 fprintf(stderr,"[external] database=%d queries=%d; graph excludes queries\n",D.n,external_nq);
}

static void usage() {
    fprintf(stderr,
        "SJG: Skyline-based Joint Graph for multi-metric nearest-neighbor search\n"
        "Usage:\n"
        "  sjg build  <database.txt> <graph.ind> [graph.sjg]\n"
        "  sjg search <database.txt> <graph.sjg> <queries.txt>\n"
        "  sjg brute  <database.txt> <queries.txt>\n"
        "Options (environment): HJG_THREADS=12, SJG_INDEX_K=15, SKY_FLOOR=16,\n"
        "  EF=128 or EF_LIST=16,32,64,128, OUT=answers.txt, OUT_DIST=1,\n"
        "  SJG_EDGE_BOUNDS=1 enables EF triangle-inequality filtering (default: off),\n"
        "  EXTERNAL_QUERY_DATA=query_objects.txt. See docs/FORMAT.md.\n");
}

static void run_queries(const Dataset& D, const Query& Q, Engine* eng, const Graph* G, const SearchOpt& opt, bool use_brute, const char* outpath) {
    std::ofstream ou(outpath);
    if (!ou) { fprintf(stderr, "[sjg] cannot write output: %s\n", outpath); exit(1); }
    for (size_t ki = 0; ki < Q.kvals.size(); ki++) {
        int k = Q.kvals[ki];
        g_dc_reset();
        double wt0 = omp_get_wtime();
        for (size_t qj = 0; qj < Q.pos.size(); qj++) {
            int raw=Q.pos[qj]; if(external_nq && (raw<0 || raw>=external_nq)) {fprintf(stderr,"invalid external query ID\n");exit(2);}
            int qp = external_nq ? D.n+raw : raw%D.n;
            Result r = use_brute ? brute(D, qp, Q.w, k)
                     : (opt.ef > 0 ? eng->search_ef(G && !external_nq ? G->perm[qp] : qp, Q.w, k, opt.ef)
                                   : eng->search(G && !external_nq ? G->perm[qp] : qp, Q.w, k, opt));
            double mx = 0;
            for (double d : r.dists) if (d > mx) mx = d;
            ou << "k " << mx << " " << r.ids.size();
            for (int id : r.ids) ou << " " << (G ? G->rid[id] : id);
            ou << "\n";
            if (getenv("OUT_DIST")) { ou << "D"; for (double d : r.dists) ou << " " << d; ou << "\n"; }
        }
        ou << "t\t" << (omp_get_wtime() - wt0) * 1000.0 / (double)Q.pos.size() << "\n";
        ou << "dc\t" << (double)g_dc_sum() / (double)Q.pos.size() << "\n";
    }
}

int main(int argc, char** argv) {
    if (argc < 3) { usage(); return 1; }
    std::string mode = argv[1];
    if (!getenv("HJG_SKYLINE")) setenv("HJG_SKYLINE", "1", 0);
    if (!getenv("SKY_FLOOR")) setenv("SKY_FLOOR", "16", 0);
    if (!getenv("HJG_THREADS")) setenv("HJG_THREADS", "1", 0);
    if (!getenv("EF") && !getenv("EF_LIST")) setenv("EF", "128", 0);
    Dataset D;
    if (!load_dataset(argv[2], D)) return 1;
    fprintf(stderr, "[sjg] data: n=%d m=%d nstr=%d\n", D.n, D.m, D.nstr);

    SearchOpt opt;
    opt.lazy_dw = !getenv("NOLAZY");
    opt.lazy_lb = getenv("LAZY_LB") != NULL;
    if (getenv("EPSI")) opt.epsi = atof(getenv("EPSI"));
    if (getenv("EF")) opt.ef = atoi(getenv("EF"));
    const char* outpath = getenv("OUT") ? getenv("OUT") : "ans-tbg.txt";

    if (getenv("NONORM")) g_use_norm = false;
    if (getenv("NOMYERS")) g_use_myers = false;
    if (mode == "testedit") {   // Myers vs DP 自测
        extern int edit_dp(const char*, int, const char*, int);
        extern int edit_myers(const char*, int, const char*, int);
        srand(7);
        int bad = 0;
        for (int t = 0; t < 50000; t++) {
            int l1 = rand() % 80, l2 = rand() % 80;
            std::string a, b;
            int alpha = (t % 2) ? 4 : 26;
            for (int i = 0; i < l1; i++) a += (char)('a' + rand() % alpha);
            for (int i = 0; i < l2; i++) b += (char)('a' + rand() % alpha);
            int d1 = edit_dp(a.c_str(), l1, b.c_str(), l2);
            int d2 = edit_myers(a.c_str(), l1, b.c_str(), l2);
            if (d1 != d2) { bad++; if (bad < 5) fprintf(stderr, "MISMATCH dp=%d myers=%d a=%s b=%s\n", d1, d2, a.c_str(), b.c_str()); }
        }
        fprintf(stderr, bad ? "[testedit] FAIL %d/50000\n" : "[testedit] PASS 50000/50000\n", bad);
        return bad ? 1 : 0;
    }
    if (mode == "build") {
        if (argc < 4) { usage(); return 1; }
        precompute_norms(D);   // cos 范数缓存: 值逐位一致(开方截断同老核), NN-descent 免 2/3 FLOPs
        std::string osjg = (argc >= 5) ? argv[4] : (std::string(argv[3]) + ".sjg");
        if (getenv("NO_NATIVE_OUT")) osjg = "";   // 高m时native边车巨大(m=50: 数百GB), 不落盘; 搜索端EF_LIST内存重建一次
        return build_graph(D, argv[3], osjg) ? 0 : 1;
    }
    if (mode == "brute") {
        if (argc < 4) { usage(); return 1; }
        Query Q;
        if (!load_query(argv[3], D.m, Q)) return 1;
        precompute_norms(D);
        external_append(D);
        run_queries(D, Q, nullptr, nullptr, opt, true, outpath);
    } else if (mode == "search") {
        if (argc < 5) { usage(); return 1; }
        if (argc >= 6) opt.epsi = atof(argv[5]);
        Graph G;
        std::string gp = argv[3];
        bool is_sjg = gp.size() > 4 && gp.substr(gp.size() - 4) == ".sjg";
        if (is_sjg) { if (!load_native(gp, D, G)) return 1; }
        else {
            const char* flag = getenv("HJG_FLATTEN_LEVELS");
            bool flatten = flag && std::string(flag) == "1";
            if (flatten && opt.ef <= 0 && !getenv("EF_LIST")) {
                fprintf(stderr, "[sjg] HJG_FLATTEN_LEVELS requires EF or EF_LIST\n");
                return 1;
            }
            if (!load_hjg_ind(gp, D, G, flatten)) return 1;
        }
        Query Q;
        if (!load_query(argv[4], D.m, Q)) return 1;
        if (!getenv("NOREORDER")) locality_reorder(D, G);
        if (getenv("NONORM")) g_use_norm = false;
        if (getenv("NOMYERS")) g_use_myers = false;
        precompute_norms(D);
        if ((getenv("NATIVE") || opt.ef > 0 || getenv("EF_LIST")) && !G.native) build_native(D, G);   // EF 模式需要 native 全维元组
        external_append(D);
        Engine eng;
        eng.init(D, G);
        if (const char* el = getenv("EF_LIST")) {   // 进程内扫多档EF(native只建一次): 每档答案写 OUT.<ef>
            char* dup = strdup(el);
            for (char* tok = strtok(dup, ","); tok; tok = strtok(nullptr, ",")) {
                SearchOpt o2 = opt; o2.ef = atoi(tok);
                std::string op = std::string(outpath) + "." + tok;
                run_queries(D, Q, &eng, &G, o2, false, op.c_str());
            }
            free(dup);
        } else
        run_queries(D, Q, &eng, &G, opt, false, outpath);
    } else { usage(); return 1; }
    return 0;
}
