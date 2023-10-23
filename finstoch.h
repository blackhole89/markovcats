/* an enumerable domain for finstoch */
#ifndef _FINSTOCH_H_
#define _FINSTOCH_H_

#include "mcat.h"
#include<map>
#include<unordered_map>
#include<functional>
#include<set>
#include<armadillo>

using namespace arma;

struct FinStochMor;
struct FinStochObj : public Obj {
    // ASSUME there are no zeroes in here
    std::vector<unsigned int> sizes;

    virtual Mor* id();
    virtual int tensor_length() {
        return sizes.size();
    }
    virtual Mor* copy_n(int n);

    virtual bool is_equal(Obj *other) {
        // assume conversion goes through
        FinStochObj *fsother = (FinStochObj*)other;

        return (sizes == fsother->sizes);
    }
    virtual FinStochObj* tensor(Obj *other) {
        // assume conversion goes through
        FinStochObj *fsother = (FinStochObj*)other;

        FinStochObj *ret = new FinStochObj();

        ret->sizes.insert( ret->sizes.end(), sizes.begin(), sizes.end() );
        ret->sizes.insert( ret->sizes.end(), fsother->sizes.begin(), fsother->sizes.end() );

        return ret;
    }
    virtual Mor* term();
    virtual Mor* select(int n);

    FinStochObj() {}
    FinStochObj(unsigned int _size) { if(_size) sizes.push_back(_size); }

    unsigned int total_dim() {
        unsigned int r=1;
        for(int i=0;i<sizes.size();++i) r*=sizes[i];
        return r;
    }

    virtual void debug_print() {
        printf("FinStochObj %lX: ", this);
        if(!sizes.size()) printf("I\n");
        else {
            printf("[%d]",sizes[0]);
            for(int i=1;i<sizes.size();++i) printf("x[%d]",sizes[i]);
            printf("\n");
        }
    }

    int tuple2coord(std::vector<int> v) {
        if(v.size() != sizes.size()) return -1;

        int stride=1;
        int ret=0;
        for(int i=v.size()-1;i>=0;--i) {
            ret += stride * v[i];
            stride *= sizes[i];
        }
        return ret;
    }
    // deterministic morphism to tuple
    FinStochMor *det_const(std::vector<int>);
};

struct FinStochMor : public Mor {
    FinStochObj *o_dom, *o_cod;
    mat data;

    virtual FinStochObj* dom() { return o_dom; }
    virtual FinStochObj* cod() { return o_cod; }
    virtual FinStochMor* cond(int n) {
        FinStochMor *ret = new FinStochMor();

        ret->o_dom = new FinStochObj(*o_dom);
        ret->o_cod = new FinStochObj(*o_cod);
        ret->o_dom->sizes.insert(ret->o_dom->sizes.end(), o_cod->sizes.end()-n, o_cod->sizes.end());
        ret->o_cod->sizes.resize(o_cod->sizes.size() - n);
        ret->data = zeros(ret->o_dom->total_dim(), ret->o_cod->total_dim());

        // in: A -> C (x) B     out: A (x) B -> C
        int adim = o_dom->total_dim();
        int bdim = ret->o_dom->total_dim()/adim;
        int cdim = ret->o_cod->total_dim();
//        printf("dims: %d %d %d\n", adim,bdim,cdim); fflush(stdout);

        for(int a=0;a<adim;++a) {
            for(int b=0;b<bdim;++b) {
                double csum = 0.0;
                for(int c=0;c<cdim;++c) {
                    csum += data(a, bdim*c+b);
                }
                if(csum==0.0) {
                    // avoid divide-by-zero when all are 0
                    for(int c=0;c<cdim;++c) {
                        ret->data(bdim*a+b, c) = 1.0/cdim;
                    }
                } else {
                    for(int c=0;c<cdim;++c) {
                        ret->data(bdim*a+b, c) = data(a, bdim*c+b)/csum;
                    }
                }
            }
        }

        return ret;
    }
    virtual FinStochMor* compose(Mor* other) {
        if(!other->dom()->is_equal(cod())) return NULL;

        // assume conversion goes through
        FinStochMor *fsother = (FinStochMor*)other;

        FinStochMor *ret = new FinStochMor();
        ret->o_dom = o_dom;
        ret->o_cod = fsother->o_cod;
        ret->data = data * fsother->data;
        return ret;
    }
    virtual FinStochMor* tensor(Mor *other) {
        FinStochMor *ret = new FinStochMor();

        // assume conversion goes through
        FinStochMor *fsother = (FinStochMor*)other;

        ret->o_dom = o_dom->tensor(fsother->o_dom);
        ret->o_cod = o_cod->tensor(fsother->o_cod);
        ret->data = kron(data, fsother->data);

        return ret;
    }

    virtual void debug_print() {
        printf("== FinStochMor %lX ==\n", this);
        printf("dom: ");
        o_dom->debug_print();
        printf("cod: ");
        o_cod->debug_print();
        data.print(std::cout);
    }
};

Mor* FinStochObj::id()
{
    FinStochMor *ret = new FinStochMor();
    ret->o_dom = ret->o_cod = this;
    ret->data = eye(total_dim(), total_dim());
    return ret;
}

Mor* FinStochObj::term()
{
    FinStochMor *ret = new FinStochMor();
    ret->o_dom = this;
    ret->o_cod = new FinStochObj();
    ret->data = ones(total_dim(), 1);
    return ret;
}

Mor* FinStochObj::copy_n(int n)
{
    FinStochMor *ret = new FinStochMor();
    ret->o_dom = this;
    ret->o_cod = this;
    int vs = total_dim();
    for(int i=1;i<n;++i) {
        ret->o_cod = ret->o_cod->tensor(this);
        vs *= total_dim();
    }
    ret->data = zeros( total_dim(), vs );

    int geosum = (vs-1)/(total_dim()-1);
    vs /= total_dim();
    for(int i=0;i<total_dim();++i) {
        ret->data(i,i*geosum)=1;
    }
    return ret;
}

Mor* FinStochObj::select(int n)
{
    // optimised, but could be constructed from term and id
    FinStochMor *ret = new FinStochMor();
    ret->o_dom = this;
    ret->o_cod = new FinStochObj(sizes[n]);

    unsigned int diml=1,dimr=1;
    for(int i=0;i<n;++i) {
        diml *= sizes[i];
    }
    for(int i=(n+1);i<sizes.size();++i) {
        dimr *= sizes[i];
    }
    ret->data = kron(kron(ones(diml,1),eye(sizes[n],sizes[n])),ones(dimr,1));
    return ret;
}

FinStochMor *FinStochObj::det_const(std::vector<int> v)
{
    FinStochMor *ret = new FinStochMor();
    ret->o_dom = new FinStochObj();
    ret->o_cod = this;
    ret->data = zeros( 1, total_dim() );
    ret->data( 0, tuple2coord(v) ) = 1;
    return ret;
}

#endif
