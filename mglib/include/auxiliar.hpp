#ifndef AUXILIAR_H_INCLUDED
#define AUXILIAR_H_INCLUDED

#include "main.hpp"
#include <time.h>

template <class T>
inline void print_vector(const vector<T>& vec)
{
    cout << "vector size: " << vec.size() << endl;
    for (typename vector<T>::const_iterator iter = vec.begin(); iter != vec.end(); iter++)
        cout << *iter << endl;
}

template <class T>
void copy_vector(vector<T>& from, vector<T>& to)
{
    to.clear();
    to.resize(from.size());
    copy(from.begin(), from.end(), to.begin());
}

template <class T>
void attach_vector(vector<T>& first, vector<T>& second)
{
    for (typename vector<T>::iterator lit = second.begin(); lit != second.end(); lit++)
        first.push_back(*lit);
}

inline void reset_time()
{
    prev_time = get_cpu_time();
}

inline double get_used_time()
{
    current_time = get_cpu_time();
    double diff = current_time - prev_time;
    prev_time = current_time;
    return diff;
}


#endif