/*
 * End-to-end parity test: the bridge (file-mem path AND shared-store path)
 * must produce the SAME overlaps as the original hifiasm CLI.
 *
 * The test is self-contained: it generates a synthetic overlapping read set,
 * writes it to a temp FASTA, runs the freshly built `hifiasm` binary on it to
 * get the ground-truth PAF (both alignment-filtered and raw-candidate modes),
 * then runs the in-memory bridge paths over the same reads and compares.
 *
 * We compare the core PAF columns (qname, q_start, q_end, strand, tname,
 * t_start, t_end, matches, block_len) as SETS, ignoring the cg:Z: CIGAR tag
 * which the in-memory sink does not carry.
 *
 * For each mode we assert:  file-mem path == CLI  AND  store path == CLI.
 *
 * Build/run (from the submodule root, after `make` so the CLI exists):
 *   make test_cli_parity && ./test_cli_parity
 */
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include <unistd.h>

#include "hifiasm_overlaps.h"

static int g_fail = 0;
#define CHECK(c,m) do{ if(!(c)){ std::fprintf(stderr,"FAIL: %s\n",(m)); ++g_fail; } }while(0)

struct Key {
    std::string q,t; uint32_t qs,qe,ts,te,nm,bl; uint8_t ss;
    bool operator<(const Key&o)const{
        if(q!=o.q)return q<o.q; if(t!=o.t)return t<o.t;
        if(qs!=o.qs)return qs<o.qs; if(qe!=o.qe)return qe<o.qe;
        if(ts!=o.ts)return ts<o.ts; if(te!=o.te)return te<o.te;
        if(nm!=o.nm)return nm<o.nm; if(bl!=o.bl)return bl<o.bl;
        return ss<o.ss;
    }
    bool operator==(const Key&o)const{
        return q==o.q&&t==o.t&&qs==o.qs&&qe==o.qe&&ts==o.ts&&te==o.te&&
               nm==o.nm&&bl==o.bl&&ss==o.ss;
    }
};

/* Parse a PAF file into the comparable key set (core columns only). */
static std::set<Key> parse_paf(const char* path){
    std::set<Key> s; std::ifstream in(path); std::string line;
    while(std::getline(in,line)){
        if(line.empty()) continue;
        std::istringstream ss(line);
        std::string qn,ql,qs,qe,st,tn,tl,ts,te,mt,bl;
        ss>>qn>>ql>>qs>>qe>>st>>tn>>tl>>ts>>te>>mt>>bl;
        Key k; k.q=qn; k.t=tn;
        k.qs=(uint32_t)strtoul(qs.c_str(),0,10); k.qe=(uint32_t)strtoul(qe.c_str(),0,10);
        k.ts=(uint32_t)strtoul(ts.c_str(),0,10); k.te=(uint32_t)strtoul(te.c_str(),0,10);
        k.nm=(uint32_t)strtoul(mt.c_str(),0,10); k.bl=(uint32_t)strtoul(bl.c_str(),0,10);
        k.ss=(st=="+")?1:0;
        s.insert(k);
    }
    return s;
}

struct Read{ std::string name,seq; };
static std::vector<Read> read_fasta(const char* path){
    std::vector<Read> v; std::ifstream in(path); std::string line; Read cur; bool have=false;
    while(std::getline(in,line)){
        if(line.empty()) continue;
        if(line[0]=='>'){ if(have)v.push_back(cur); cur=Read(); cur.name=line.substr(1); have=true; }
        else cur.seq+=line;
    }
    if(have)v.push_back(cur);
    return v;
}

static std::string name_of(const char*n,const uint64_t*o,uint64_t i){ return std::string(n+o[i],n+o[i+1]); }

static std::set<Key> to_keys(const hifiasm_overlap_t*ov,uint64_t n,const char*nm,const uint64_t*off){
    std::set<Key> s;
    for(uint64_t i=0;i<n;i++){
        Key k; k.q=name_of(nm,off,ov[i].q_id); k.t=name_of(nm,off,ov[i].t_id);
        k.qs=ov[i].q_start; k.qe=ov[i].q_end; k.ts=ov[i].t_start; k.te=ov[i].t_end;
        k.nm=ov[i].n_match; k.bl=ov[i].block_len; k.ss=ov[i].is_same_strand;
        s.insert(k);
    }
    return s;
}

static void run_and_check(const char* fasta, const char* cli_paf,
                          int raw_candidates, const char* tag){
    std::set<Key> oracle = parse_paf(cli_paf);
    std::fprintf(stderr,"[%s] CLI overlaps=%zu\n", tag, oracle.size());
    CHECK(!oracle.empty(), "CLI produced overlaps");

    std::vector<Read> reads = read_fasta(fasta);

    hifiasm_ovlp_opt_t opt; std::memset(&opt,0,sizeof(opt));
    opt.threads=4; opt.raw_candidates=raw_candidates;

    /* ---- file-mem path ---- */
    const char* files[1]={fasta};
    hifiasm_overlap_t* fov=0; uint64_t fn=0; char* fnm=0; uint64_t* foff=0; uint64_t fnr=0;
    int rc=hifiasm_detect_overlaps_mem(files,1,&opt,&fov,&fn,&fnm,&foff,&fnr,NULL,NULL);
    CHECK(rc==0,"file-mem rc==0");
    std::set<Key> fset=to_keys(fov,fn,fnm,foff);
    std::fprintf(stderr,"[%s] file-mem overlaps=%zu\n", tag, fset.size());
    CHECK(fset==oracle, "file-mem path == CLI");

    /* ---- store path ---- */
    std::vector<hifiasm_read_t> hr(reads.size());
    for(size_t i=0;i<reads.size();i++){
        hr[i].seq=reads[i].seq.data(); hr[i].seq_len=reads[i].seq.size();
        hr[i].name=reads[i].name.data(); hr[i].name_len=(uint32_t)reads[i].name.size();
    }
    rc=hifiasm_reads_store_load(hr.data(),hr.size());
    CHECK(rc==0,"store load rc==0");
    hifiasm_overlap_t* sov=0; uint64_t sn=0; char* snm=0; uint64_t* soff=0; uint64_t snr=0;
    rc=hifiasm_detect_overlaps_from_store(&opt,&sov,&sn,&snm,&soff,&snr,NULL,NULL);
    CHECK(rc==0,"store rc==0");
    std::set<Key> sset=to_keys(sov,sn,snm,soff);
    std::fprintf(stderr,"[%s] store overlaps=%zu\n", tag, sset.size());
    CHECK(sset==oracle, "store path == CLI");
    CHECK(snr==reads.size(), "store read count matches input");

    hifiasm_overlaps_mem_free(fov,fnm,foff,NULL);
    hifiasm_overlaps_mem_free(sov,snm,soff,NULL);
    hifiasm_reads_store_release();
}

/* Reproducible pseudo-random ACGT reference (xorshift32). */
static std::string random_ref(size_t n,uint32_t seed){
    static const char b[4]={'A','C','G','T'}; std::string s(n,'A'); uint32_t x=seed?seed:1;
    for(size_t i=0;i<n;i++){ x^=x<<13; x^=x>>17; x^=x<<5; s[i]=b[x&3]; }
    return s;
}

/* Cut the reference into overlapping windows with a few substitution errors,
 * so real overlaps exist and the alignment filter has work to do. */
static std::vector<Read> gen_reads(){
    std::vector<Read> v; std::string ref=random_ref(60000,42);
    const size_t win=8000, step=2500; uint32_t x=99; int idx=0;
    for(size_t s=0;s+win<=ref.size();s+=step){
        std::string seq=ref.substr(s,win);
        for(int e=0;e<20;e++){ x^=x<<13; x^=x>>17; x^=x<<5;
            size_t p=x%seq.size(); seq[p]="ACGT"[(x>>8)&3]; }
        Read r; r.name="read"+std::to_string(idx++); r.seq=seq; v.push_back(r);
    }
    /* ambiguous bases in one read to exercise the N_site path */
    if(v.size()>1) for(size_t i=200;i<210&&i<v[1].seq.size();i++) v[1].seq[i]='N';
    return v;
}

static bool write_fasta(const std::vector<Read>&v,const char*p){
    FILE*f=std::fopen(p,"w"); if(!f)return false;
    for(const auto&r:v) std::fprintf(f,">%s\n%s\n",r.name.c_str(),r.seq.c_str());
    std::fclose(f); return true;
}

int main(int argc,char**argv){
    /* Locate the hifiasm CLI (built by `make`); allow override via argv[1]. */
    const char* cli = (argc>1)? argv[1] : "./hifiasm";
    if(access(cli,X_OK)!=0){
        std::fprintf(stderr,"SKIP: hifiasm CLI not found at %s (run `make` first)\n",cli);
        return 77; /* automake-style skip code */
    }

    std::vector<Read> reads=gen_reads();
    CHECK(reads.size()>=4,"generated several reads");

    char fa[]="/tmp/hifiasm_cli_parity_XXXXXX";
    int fd=mkstemp(fa); CHECK(fd>=0,"mkstemp"); if(fd>=0) close(fd);
    CHECK(write_fasta(reads,fa),"write fasta");

    std::string pfx=std::string(fa)+".out";
    std::string def_paf=pfx+".ovlp.paf";
    std::string raw_paf=pfx+".candidates.paf";

    char cmd[1024];
    std::snprintf(cmd,sizeof(cmd),"%s -t4 -o %s %s >/dev/null 2>&1",cli,pfx.c_str(),fa);
    CHECK(system(cmd)==0,"CLI default run");
    std::snprintf(cmd,sizeof(cmd),"%s -t4 --dbg-ovec -o %s %s >/dev/null 2>&1",cli,pfx.c_str(),fa);
    CHECK(system(cmd)==0,"CLI raw run");

    run_and_check(fa, def_paf.c_str(), /*raw*/0, "default");
    run_and_check(fa, raw_paf.c_str(), /*raw*/1, "raw");

    std::remove(fa); std::remove(def_paf.c_str()); std::remove(raw_paf.c_str());

    if(g_fail){ std::fprintf(stderr,"%d CHECK(s) FAILED\n",g_fail); return 1; }
    std::fprintf(stderr,"ALL CLI-PARITY CHECKS PASSED\n");
    return 0;
}
