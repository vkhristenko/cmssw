#ifndef testing_cmp_cpu_vs_gpu_common_h
#define testing_cmp_cpu_vs_gpu_common_h

template<int N, int S, int E>
struct comparator {
    void create(std::string const& name, edm::Service<TFileService>& fs) {
        name_ = name;
        correlation = fs->make<TH2F>(
            (name_ + "_comparison").c_str(),
            (name_ + " comparison").c_str(),
            N, S, E,
            N, S, E);
        ratio = fs->make<TH1F>(
            (name_ + "_ratio").c_str(),
            (name_ + " ratio").c_str(),
            100, -10, 10);
        v1 = fs->make<TH1F>(
            (name_ + "_cpu").c_str(),
            (name_ + " cpu").c_str(),
            N, S, E);
        v2 = fs->make<TH1F>(
            (name_ + "_gpu").c_str(),
            (name_ + " gpu").c_str(),
            N, S, E);
    }
    
    std::string name_;
    TH2F *correlation;
    TH1F *ratio;
    TH1F *v1, *v2;
};

#endif
