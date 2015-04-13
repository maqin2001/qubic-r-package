#include <iostream>
#include <fstream>

#include <queue>
#include <vector>
#include <cassert>

#include "struct.h"
#include "edge_list.h"

/*we can reduce the HEAP_SIZE when the data contain so many genes so that memory is not enough*/
#define HEAP_SIZE 20000000

int str_intersect_r(const std::vector<discrete> &s1, const std::vector<discrete> &s2) {
  assert(s1.size() == s2.size());
  int common_cnt = 0;
  /* s1 and s2 of equal length, so we check s1 only */
  for (size_t i = 0; i < s1.size(); i++)
    if ((s1[i] != 0) && (s1[i] == s2[i])) // Changed order by ZHANG Yu 
      common_cnt++;
  return common_cnt;
}

static int edge_cmpr(void * a, void * b)
{
  int score_a, score_b;
  score_a = ((Edge *)a)->score;
  score_b = ((Edge *)b)->score;

  if (score_a < score_b) return -1;
  if (score_a == score_b) return 0;
  return 1;
}

static void fh_insert_fixed(fibheap *h, Edge *data, const Edge *cur_min)
{
  if (h->fh_n < HEAP_SIZE) {
    fh_insert(h, data);
  }
  else {
    if (edge_cmpr(&cur_min, data) < 0) {
      /* Remove least value and renew */
      fh_extractmin(h);
      fh_insert(h, data);
      /* Keep a memory of the current min */
      cur_min = (Edge *)fh_min(h);
    }
  }
}

static void fh_dump(fibheap *h, std::vector<Edge *> &data_array) {
  int i;
  int n = h->fh_n;
  for (i = n - 1; i >= 0; i--)
    data_array[i] = (Edge *)fh_extractmin(h);
}

int EdgeList::get_key(const Edge* s) { return s->score - col_width; }

class MyIterator : public std::iterator < std::bidirectional_iterator_tag, const discrete >
{
  const std::vector<std::vector<discrete>> &p;
  size_t x;
  size_t y;
public:
  MyIterator(const std::vector<std::vector<discrete>> &p, size_t x = 0, size_t y = 0) :p(p), x(x), y(y) {}
  MyIterator(const MyIterator& mit) : p(mit.p), x(mit.x), y(mit.y) {}
  MyIterator& operator++() { if (x == p.size()) y++; x++; return *this; }
  MyIterator operator++(int) { MyIterator tmp(*this); operator++(); return tmp; }
  bool operator==(const MyIterator& rhs) const { return p == rhs.p; }
  bool operator!=(const MyIterator& rhs) const { return p != rhs.p; }
  const discrete& operator*() const { return p[x][y]; }
};

int get_key1(const Edge& s) { return s.score; }
class MyIterator1 : public std::iterator<std::bidirectional_iterator_tag, const Edge> {
  const std::vector<std::vector<discrete>> &p;
  size_t x;
  size_t y;
  size_t size;
public:
  MyIterator1(const std::vector<std::vector<discrete>> &p, size_t x = 0, size_t y = 1) :p(p), x(x), y(y), size(p[0].size()*(p[0].size()+1)/2) {}
  MyIterator1(const MyIterator1& mit) : MyIterator1(mit.p, mit.x, mit.y) {}
  MyIterator1& operator++() { if (y == p.size() - 1) { x++; y = x + 1; } else y++; return *this; }
  MyIterator1& operator--() { if (y == x + 1){ x--; y = p.size() - 1; } else y--; return *this; }
  MyIterator1 operator++(int) { MyIterator1 tmp(*this); operator++(); return tmp; }
  MyIterator1 operator--(int) { MyIterator1 tmp(*this); operator--(); return tmp; }
  bool operator==(const MyIterator1& rhs) const { return (p == rhs.p) && (x == rhs.x) && (y == rhs.y); }
  bool operator!=(const MyIterator1& rhs) const { return (p != rhs.p) || (x != rhs.x) || (y != rhs.y); }
  Edge& operator*() const {
    Edge edge;
    edge.gene_one = x;
    edge.gene_two = y;
    edge.score = p[0].size() - str_intersect_r(p[x], p[y]) - 1;
    return edge;
  }
  MyIterator1 begin() const { return MyIterator1(p, 0, 1); }
  MyIterator1 end() const { return MyIterator1(p, p.size() - 1, p.size()); }
};

struct CompEventByPtr {
  bool operator()(const Edge* pEvent1, const Edge* pEvent2) const {
    return pEvent1->score >= pEvent2->score;
  }
};

const std::vector<Edge *>& EdgeList::get_edge_list() const { return edge_list; }
EdgeList::EdgeList(const std::vector<std::vector<discrete>>& arr_c, int& COL_WIDTH) {
  if (COL_WIDTH == 2) COL_WIDTH = MAX(arr_c[0].size() / 20, 2);
#if 1
  Edge * edge;
  int cnt;
  int rec_num = 0;

  /* Allocating heap structure */
  fibheap *heap = fh_makeheap();
  fh_setcmp(heap, edge_cmpr);

  /* Generating seed list and push into heap */
  progress("Generating seed list (minimum weight %d)", COL_WIDTH);
  Edge __cur_min = { 0, 0, COL_WIDTH };
  Edge *_cur_min = &__cur_min;
  Edge **cur_min = &_cur_min;
  /* iterate over all genes to retrieve all edges */
  for (size_t i = 0; i < arr_c.size(); i++)
    for (size_t j = i + 1; j < arr_c.size(); j++) {
      cnt = str_intersect_r(arr_c[i], arr_c[j]);
      if (cnt < _cur_min->score) continue;
      edge = new Edge();
      edge->gene_one = i;
      edge->gene_two = j;
      edge->score = cnt;
      fh_insert_fixed(heap, edge, *cur_min);
    }
  rec_num = heap->fh_n;
  if (rec_num == 0)
    errAbort("Not enough overlap between genes");
  /* sort the seeds */
  printf("%d seeds generated\n", rec_num);
  edge_list.resize(rec_num);
  fh_dump(heap, edge_list);
  printf("%d seeds dumped\n", edge_list.size());
#else
  Edge * edge;
  int cnt;
  int rec_num = 0;

  /* Allocating heap structure */
  std::priority_queue<Edge *, std::vector<Edge *>, CompEventByPtr> q;

  /* Generating seed list and push into heap */
  progress("Generating seed list (minimum weight %d)", COL_WIDTH);
  Edge __cur_min = { 0, 0, COL_WIDTH };
  Edge *_cur_min = &__cur_min;
  Edge **cur_min = &_cur_min;
  /* iterate over all genes to retrieve all edges */
  for (size_t i = 0; i < arr_c.size(); i++)
    for (size_t j = i + 1; j < arr_c.size(); j++) {
      cnt = str_intersect_r(arr_c[i], arr_c[j]);
      if (cnt < _cur_min->score) continue;
      edge = new Edge();
      edge->gene_one = i;
      edge->gene_two = j;
      edge->score = cnt;

      if (q.size() < HEAP_SIZE) {
        q.push(edge);
      }
      else {
        if (_cur_min->score <= edge->score) {
          /* Remove least value and renew */
          q.pop();
          q.push(edge);
          /* Keep a memory of the current min */
          _cur_min = q.top();
        }
      }
    }
  rec_num = q.size();
  if (rec_num == 0)
    errAbort("Not enough overlap between genes");
  /* sort the seeds */
  printf("%d seeds generated\n", rec_num);
  edge_list.resize(rec_num);

  for (int i = rec_num - 1; i >= 0; i--) {
    edge_list[i] = q.top();
    q.pop();
  }
#endif
  //std::ofstream myfile;
  //myfile.open("seeds1.txt");
  //for (size_t i = 0; i < edge_list.size(); i++) {
  //  myfile << edge_list[i]->gene_one << " " << edge_list[i]->gene_two << " " << edge_list[i]->score << std::endl;
  //}
  //myfile.close();
}

EdgeList::~EdgeList() {
  for (size_t i = 0; i < edge_list.size(); i++)
    delete(edge_list[i]);
}