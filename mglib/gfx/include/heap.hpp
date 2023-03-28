#ifndef HHEAP_INCLUDED
#define HHEAP_INCLUDED

#include <vector>
#include <stddef.h>

// source grabbed from MixTex

class MxHeapable;
typedef std::vector<MxHeapable *> MxHeapBase;

class MxHeapable
{
private:
  float import;
  int token;

public:
  MxHeapable() { not_in_heap(); heap_key(0.0f); }

  inline bool is_in_heap() { return token != -47; }
  inline void not_in_heap() { token = -47; }
  inline int get_heap_pos() { return token; }
  inline void set_heap_pos(int t) { token=t; }

  inline void  heap_key(float k) { import=k; }
  inline float heap_key() const  { return import; }
};



class MxHeap : private MxHeapBase
{
private:
  void place(MxHeapable *x, unsigned int i);
  void swap(unsigned int i, unsigned int j);

  unsigned int parent(unsigned int i) { return (i-1)/2; }
  unsigned int left(unsigned int i) { return 2*i+1; }
  unsigned int right(unsigned int i) { return 2*i+2; }

  void upheap(unsigned int i);
  void downheap(unsigned int i);

public:
  MxHeap() { }
  MxHeap(unsigned int n) { }

  void insert(MxHeapable *t) { insert(t, t->heap_key()); }
  void insert(MxHeapable *, float);
  void update(MxHeapable *t) { update(t, t->heap_key()); }
  void update(MxHeapable *, float);
  void reset() { MxHeapBase::clear(); }

  size_t size() const { return MxHeapBase::size(); }
  MxHeapable       *item(int i)       { return (*this)[i]; }
  const MxHeapable *item(int i) const { return (*this)[i]; }
  MxHeapable *extract();
  MxHeapable *top() { return (size()<1 ? (MxHeapable *)NULL : item(0)); }
  MxHeapable *remove(MxHeapable *);
};

#endif
