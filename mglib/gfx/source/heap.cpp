/************************************************************************

  source grabbed from MixTex

 ************************************************************************/

#include "heap.hpp"

#define ref(i) ((*this)[i])


inline void MxHeap::place(MxHeapable *x, unsigned int i)
{
  ref(i) = x;
  x->set_heap_pos(i);
}

void MxHeap::swap(unsigned int i, unsigned int j)
{
  MxHeapable *tmp = ref(i);

  place(ref(j), i);
  place(tmp, j);
}

void MxHeap::upheap(unsigned int i)
{
  MxHeapable *moving = ref(i);
  unsigned int index = i;
  unsigned int p = parent(i);

  while( index>0 && moving->heap_key() > ref(p)->heap_key() )
  {
    place(ref(p), index);
    index = p;
    p = parent(p);
  }

  if( index != i )
    place(moving, index);
}

void MxHeap::downheap(unsigned int i)
{
  MxHeapable *moving = ref(i);
  unsigned int index = i;
  unsigned int l = left(i);
  unsigned int r = right(i);
  unsigned int largest;

  while( l<size() )
  {
    if( r<size() && ref(l)->heap_key() < ref(r)->heap_key() )
      largest = r;
    else
      largest = l;

    if( moving->heap_key() < ref(largest)->heap_key() )
    {
      place(ref(largest), index);
      index = largest;
      l = left(index);
      r = right(index);
    }
    else
      break;
  }

  if( index != i )
    place(moving, index);
}

////////////////////////////////////////////////////////////////////////
//
// Exported interface to the heap
//

void MxHeap::insert(MxHeapable *t, float v)
{
  t->heap_key(v);

  push_back(t);
  unsigned int i = size()-1;
  t->set_heap_pos(i);

  upheap(i);
}

void MxHeap::update(MxHeapable *t, float v)
{
  t->heap_key(v);

  unsigned int i = t->get_heap_pos();

  if( i>0 && v>ref(parent(i))->heap_key() )
    upheap(i);
  else
    downheap(i);
}

MxHeapable *MxHeap::extract()
{
  if( size() < 1 ) return NULL;

  swap(0, size()-1);
  MxHeapable *dead = ref(size()-1);
  pop_back();

  downheap(0);
  dead->not_in_heap();
  return dead;
}

MxHeapable *MxHeap::remove(MxHeapable *t)
{
  if( !t->is_in_heap() ) return NULL;

  int i = t->get_heap_pos();
  swap(i, size()-1);
  pop_back();
  t->not_in_heap();

  if( ref(i)->heap_key() < t->heap_key() )
    downheap(i);
  else
    upheap(i);

  return t;
}
