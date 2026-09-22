#include <flint/nmod_mat.h>
#include <stdint.h>

static void read_mat(nmod_mat_t M, const uint64_t *a, long r, long c, unsigned long q) {
    nmod_mat_init(M, r, c, q);
    for (long i=0;i<r;i++) for (long j=0;j<c;j++) nmod_mat_entry(M,i,j)=a[i*c+j]%q;
}
long ff_rank(const uint64_t *a,long r,long c,unsigned long q) {
    nmod_mat_t M; read_mat(M,a,r,c,q);
    long k=nmod_mat_rank(M); nmod_mat_clear(M); return k;
}
long ff_rref(const uint64_t *a,long r,long c,unsigned long q,uint64_t *out) {
    nmod_mat_t M; read_mat(M,a,r,c,q);
    long k=nmod_mat_rref(M);
    for(long i=0;i<r;i++) for(long j=0;j<c;j++) out[i*c+j]=nmod_mat_entry(M,i,j);
    nmod_mat_clear(M); return k;
}
long ff_kernel(const uint64_t *a,long r,long c,unsigned long q,uint64_t *out) {
    nmod_mat_t M,K; read_mat(M,a,r,c,q); nmod_mat_init(K,c,c,q);
    long k=nmod_mat_nullspace(K,M);
    for(long i=0;i<k;i++) for(long j=0;j<c;j++) out[i*c+j]=nmod_mat_entry(K,j,i);
    nmod_mat_clear(K); nmod_mat_clear(M); return k;
}
uint64_t ff_det(const uint64_t *a,long r,unsigned long q) {
    nmod_mat_t M; read_mat(M,a,r,r,q);
    uint64_t d=nmod_mat_det(M); nmod_mat_clear(M); return d;
}
int ff_solve(const uint64_t *a,long r,long c,const uint64_t *b,long d,unsigned long q,uint64_t *out) {
    nmod_mat_t A,B,X;
    read_mat(A,a,r,c,q); read_mat(B,b,r,d,q); nmod_mat_init(X,c,d,q);
    int ok=nmod_mat_can_solve(X,A,B);
    if(ok) for(long i=0;i<c;i++) for(long j=0;j<d;j++) out[i*d+j]=nmod_mat_entry(X,i,j);
    nmod_mat_clear(X); nmod_mat_clear(B); nmod_mat_clear(A); return ok;
}
void ff_mul(const uint64_t *a,long r,long c,const uint64_t *b,long d,unsigned long q,uint64_t *out) {
    nmod_mat_t A,B,C;
    read_mat(A,a,r,c,q); read_mat(B,b,c,d,q); nmod_mat_init(C,r,d,q);
    nmod_mat_mul(C,A,B);
    for(long i=0;i<r;i++) for(long j=0;j<d;j++) out[i*d+j]=nmod_mat_entry(C,i,j);
    nmod_mat_clear(C); nmod_mat_clear(B); nmod_mat_clear(A);
}
