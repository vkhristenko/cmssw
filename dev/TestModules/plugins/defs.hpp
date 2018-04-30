#ifndef dev_defs_h
#define dev_defs_h

#include <string>
#include <iostream>
#include <mutex>
#include <sstream>

/*
namespace par {

struct cout : public std::stringstream {
    static inline std::mutex cout_mutex;
    ~cout() {
        std::lock_guard<std::mutex> l{cout_mutex};
        std::cout << rdbuf();
    }
};

}*/

inline void print(std::string const& str) {
    std::cout << str << "\n";
}

#endif
