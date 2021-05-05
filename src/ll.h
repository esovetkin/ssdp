#ifndef ll_h
#define ll_h
#include <stdlib.h>
#include <stdbool.h>

typedef struct ll_node ll_node;
struct ll_node {
	void *data;
	ll_node *prev;
	ll_node *next;
};

#define DATA(n_n_n) ((n_n_n)->data)
#define NEXT(n_n_n) ((n_n_n)->next)
#define PREV(n_n_n) ((n_n_n)->prev)

ll_node *ll_first(ll_node *n);
ll_node *ll_last(ll_node *n);
size_t ll_length(const ll_node *n);
ll_node *ll_insert_before(ll_node *n, void *data);
ll_node *ll_insert_after(ll_node *n, void *data);
ll_node *ll_cut_before(ll_node *n);
ll_node *ll_cut_after(ll_node *n);
void ll_glue(ll_node *prev, ll_node *next);
void *ll_remove(ll_node *n);
ll_node *ll_find(ll_node *n, bool (*predicate)(void *));
ll_node *ll_rfind(ll_node *n, bool (*predicate)(void *));
ll_node *ll_find_r(ll_node *n, bool (*predicate)(void *, void *), void *extra);
ll_node *ll_rfind_r(ll_node *n, bool (*predicate)(void *, void *), void *extra);
void ll_map(ll_node *n, void *(*func) (void*));
void ll_map_r(ll_node *n, void *(*func) (void*, void*), void *extra);
void ll_nodemap(ll_node *n, ll_node *(*func) (ll_node *));
void ll_nodemap_r(ll_node *n, ll_node *(*func) (ll_node *, void *), void *extra);
void ll_destroy(ll_node *n);
void ll_mapdestroy(ll_node *n, void (func) (void*));
void ll_mapdestroy_r(ll_node *n, void (func) (void*, void*), void *state);

#ifdef LL_CIRCULAR
bool ll_is_circular(const ll_node *n);
size_t ll_clength(const ll_node *n);
ll_node *ll_cfind(ll_node *n, bool (*predicate)(void *));
ll_node *ll_crfind(ll_node *n, bool (*predicate)(void *));
ll_node *ll_cfind_r(ll_node *n, bool (*predicate)(void *, void *), void *extra);
ll_node *ll_crfind_r(ll_node *n, bool (*predicate)(void *, void *), void *extra);
void ll_cmap(ll_node *n, void *(*func) (void*));
void ll_cmap_r(ll_node *n, void *(*func) (void*, void*), void *extra);
#endif

#endif
