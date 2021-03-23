#include "ll.h"
ll_node *ll_first(ll_node *n) { /*{{{*/
	if (n == NULL) { return NULL; }
	for (; n->prev!=NULL; n=n->prev);
	return n;
} /*}}}*/
ll_node *ll_last(ll_node *n) { /*{{{*/
	if (n == NULL) { return NULL; }
	for (; n->next!=NULL; n=n->next);
	return n;
} /*}}}*/
size_t ll_length(const ll_node *n) { /*{{{*/
	if (n == NULL) { return 0; }
	size_t c = 0;
	while (n != NULL) {
		++c;
		n = n->next;
	}
	return c;
} /*}}}*/
ll_node *ll_insert_before(ll_node *n, void *data) { /*{{{*/
	ll_node *o = malloc(sizeof(ll_node));
	if (o == NULL) { return NULL; }
	o->data = data;
	if (n == NULL) { /*{{{*/
		o->prev = NULL;
		o->next = NULL;
	} /*}}}*/
	else { /*{{{*/
		if (n->prev != NULL) {
			n->prev->next = o;
		}
		o->prev = n->prev;
		o->next = n;
		n->prev = o;
	} /*}}}*/
	return o;
} /*}}}*/
ll_node *ll_insert_after(ll_node *n, void *data) { /*{{{*/
	ll_node *o = malloc(sizeof(ll_node));
	if (o == NULL) { return NULL; }
	o->data = data;
	if (n == NULL) { /*{{{*/
		o->prev = NULL;
		o->next = NULL;
	} /*}}}*/
	else { /*{{{*/
		if (n->next != NULL) {
			n->next->prev = o;
		}
		o->prev = n;
		o->next = n->next;
		n->next = o;
		if (n->prev == n) {
			n->prev = o;
		}
	} /*}}}*/
	return o;
} /*}}}*/
ll_node *ll_cut_before(ll_node *n) { /*{{{*/
	if (n == NULL) {
		return NULL;
	}
	ll_node *prev = n->prev;
	if (prev==NULL) {
		return NULL;
	}
	else {
		n->prev = NULL;
		prev->next = NULL;
		return prev;
	}
} /*}}}*/
ll_node *ll_cut_after(ll_node *n) { /*{{{*/
	if (n == NULL) {
		return NULL;
	}
	ll_node *next = n->next;
	if (next==NULL) {
		return NULL;
	}
	else {
		n->next = NULL;
		next->prev = NULL;
		return next;
	}
} /*}}}*/
void ll_glue(ll_node *prev, ll_node *next) { /*{{{*/
	if (prev != NULL) {
		prev->next = next;
	}
	if(next != NULL) {
		next->prev = prev;
	}
} /*}}}*/
ll_node *ll_find(ll_node *n, bool (*predicate)(void *)) { /*{{{*/
	for (; n!=NULL; n=n->next) {
		if (predicate(n->data)) {
			return n;
		}
	}
	return NULL;
} /*}}}*/
ll_node *ll_rfind(ll_node *n, bool (*predicate)(void *)) { /*{{{*/
	for (; n!=NULL; n=n->prev) {
		if (predicate(n->data)) {
			return n;
		}
	}
	return NULL;
} /*}}}*/
ll_node *ll_find_r(ll_node *n, bool (*predicate)(void *, void *), void *extra) { /*{{{*/
	for (; n!=NULL; n=n->next) {
		if (predicate(n->data, extra)) {
			return n;
		}
	}
	return NULL;
} /*}}}*/
ll_node *ll_rfind_r(ll_node *n, bool (*predicate)(void *, void *), void *extra) { /*{{{*/
	for (; n!=NULL; n=n->prev) {
		if (predicate(n->data, extra)) {
			return n;
		}
	}
	return NULL;
} /*}}}*/
void *ll_remove(ll_node *n) { /*{{{*/
	if (n == NULL) { return NULL; }
	if (n->prev != NULL) {
		n->prev->next = n->next;
	}
	if (n->next != NULL) {
		n->next->prev = n->prev;
	}
	void *data = n->data;
	free(n);
	return data;
} /*}}}*/
void ll_map(ll_node *n, void *(*func) (void*)) { /*{{{*/
	if (n == NULL) { return; }
	for (; n!=NULL; n=n->next) {
		n->data = func(n->data);
	}
} /*}}}*/
void ll_map_r(ll_node *n, void *(*func) (void*, void*), void *extra) { /*{{{*/
	for (; n!=NULL; n=n->next) {
		n->data = func(n->data, extra);
	}
} /*}}}*/
void ll_nodemap(ll_node *n, ll_node *(*func) (ll_node *)) { /*{{{*/
	while (n != NULL) {
		n = func(n);
	}
} /*}}}*/
void ll_nodemap_r(ll_node *n, ll_node *(*func) (ll_node *, void *), void *extra) { /*{{{*/
	while (n != NULL) {
		n = func(n, extra);
	}
} /*}}}*/
void ll_destroy(ll_node *n) { /*{{{*/
	ll_cut_before(n);
	while (n != NULL) {
		ll_node *o = n;
		n = n->next;
		free(o);
	}
} /*}}}*/
void ll_mapdestroy(ll_node *n, void (func) (void*)) { /*{{{*/
	ll_cut_before(n);
	while (n != NULL) {
		ll_node *o = n;
		func(n->data);
		n = n->next;
		free(o);
	}
} /*}}}*/
void ll_mapdestroy_r(ll_node *n, void (func) (void*, void*), void *state) { /*{{{*/
	ll_cut_before(n);
	while (n != NULL) {
		ll_node *o = n;
		func(n->data, state);
		n = n->next;
		free(o);
	}
} /*}}}*/
#ifdef LL_CIRCULAR
bool ll_is_circular(const ll_node *n) { /*{{{*/
	const ll_node *o = n;
	while (true) {
		if (n == NULL) {
			return false;
		}
		n = n->next;
		if (n == o) {
			return true;
		}
	}
}/*}}}*/
size_t ll_clength(const ll_node *n) { /*{{{*/
	if (n == NULL) { return 0; }
	size_t c = 0;
	const ll_node *o = n;
	do {
		++c;
		n = n->next;
	} while (n != o && n != NULL);
	return c;
} /*}}}*/
ll_node *ll_cfind(ll_node *n, bool (*predicate)(void *)) { /*{{{*/
	if (n == NULL) { return NULL; }
	ll_node *o = n;
	do {
		if (predicate(n->data)) {
			return n;
		}
		n = n->next;
	} while (n != o && n != NULL);
	return NULL;
} /*}}}*/
ll_node *ll_crfind(ll_node *n, bool (*predicate)(void *)) { /*{{{*/
	if (n == NULL) { return NULL; }
	ll_node *o = n;
	do {
		if (predicate(n->data)) {
			return n;
		}
		n = n->prev;
	} while (n != o && n != NULL);
	return NULL;
} /*}}}*/
ll_node *ll_cfind_r(ll_node *n, bool (*predicate)(void *, void *), void *extra) { /*{{{*/
	if (n == NULL) { return NULL; }
	ll_node *o = n;
	do {
		if (predicate(n->data, extra)) {
			return n;
		}
		n = n->next;
	} while (n != o && n != NULL);
	return NULL;
} /*}}}*/
ll_node *ll_crfind_r(ll_node *n, bool (*predicate)(void *, void *), void *extra) { /*{{{*/
	if (n == NULL) { return NULL; }
	ll_node *o = n;
	do {
		if (predicate(n->data, extra)) {
			return n;
		}
		n = n->prev;
	} while (n != o && n != NULL);
	return NULL;
} /*}}}*/
void ll_cmap(ll_node *n, void *(*func) (void*)) { /*{{{*/
	if (n == NULL) { return; }
	ll_node *o = n;
	do {
		n->data = func(n->data);
		n = n->next;
	} while (n != o && n != NULL);
} /*}}}*/
void ll_cmap_r(ll_node *n, void *(*func) (void*, void*), void *extra) { /*{{{*/
	if (n == NULL) { return; }
	ll_node *o = n;
	do {
		n->data = func(n->data, extra);
		n = n->next;
	} while (n != o && n != NULL);
} /*}}}*/
#endif
