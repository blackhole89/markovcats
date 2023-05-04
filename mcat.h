#include<vector>
#include<map>
#include<unordered_map>
#include<functional>
#include<set>

/* abstract interface to markov categories */
template<class M> struct Obj {
    virtual M* id()=0;
    virtual M* samp()=0;
};
template<class M, class O> struct Mor {
    virtual O* dom()=0;
    virtual O* cod()=0;
    virtual M* sharp()=0;
    virtual M* cond(std::vector<O*> os)=0;
};

template<class M, class O> struct MCat {
    virtual M* arrange(std::vector<int>, std::vector<M*>)=0;
    virtual M* compose(M*, M*)=0;
};



/* an enumerable domain for finstoch */
template<class P, int n> struct FinStochVal {
    int order;
    const int max = n;

    std::unordered_map<FinStochVal<P,n>*, P> tab;
};

template<class P, class V> struct FinStochMor;

template<class P, class V> struct FinStochObj : public Obj<FinStochMor<P,V> > {
    using Mor = FinStochMor<P,V>;
    using Obj = FinStochObj<P,V>;
    virtual Mor *id();
    virtual Mor *samp();

    std::vector<int> type; 
};

template<class P, class V> struct FinStochMor : public Mor<FinStochMor<P,V>, FinStochObj<P,V> > {
    using Mor = FinStochMor<P,V>;
    using Obj = FinStochObj<P,V>;
    virtual Obj* dom() { return src; };
    virtual Obj* cod() { return tgt; };
    virtual Mor* sharp();
    virtual Mor* cond(std::set<int> os);

    Obj *src,*tgt;
    virtual std::vector<V*> eval_on(std::vector<V*> v)=0;
};

template<class P, class V> struct FinStoch : public MCat<FinStochMor<P,V>, FinStochObj<P,V> > {
    using Mor = FinStochMor<P,V>;
    using Obj = FinStochObj<P,V>;
    virtual Mor* arrange(std::vector<int>, std::vector<Mor*>);
    virtual Mor* compose(Mor*, Mor*);
};

template<class P, class V> FinStochMor<P,V> *FinStoch<P,V>::arrange(std::vector<int>, std::vector<FinStochMor<P,V>*>){
    return NULL;
}

/* composition of morphisms */
template<class P, class V> struct FinStochMorComposition : public FinStochMor<P,V> {
    using Mor = FinStochMor<P,V>;
    Mor *m1, *m2;
    virtual std::vector<V*> eval_on(std::vector<V*> v) {
        m2->eval_on(m1->eval_on(v));
    }
};

template<class P, class V> FinStochMor<P,V> *FinStoch<P,V>::compose(FinStochMor<P,V>* m1, FinStochMor<P,V>* m2) {
    if (m1->cod()->type != m2->dom()->type) return NULL;

    FinStochMorComposition<P,V> *ret = new FinStochMorComposition<P,V>;
    ret->src=m1->src;
    ret->tgt=m2->tgt;
    ret->m1=m1;
    ret->m2=m2;
    return ret;
}


/* id and samp morphisms on an object */
template<class P, class V> struct FinStochMorId : public FinStochMor<P,V> {
    virtual std::vector<V*> eval_on(std::vector<V*> v) {
        return v;
    }
};

template<class P, class V> FinStochMor<P,V> *FinStochObj<P,V>::id() {
    FinStochMorId<P,V> *ret = new FinStochMorId<P,V>;
    ret->src=ret->tgt=this;
    return ret;
}

template<class P, class V> struct FinStochMorSamp : public FinStochMor<P,V> {
    virtual std::vector<V*> eval_on(std::vector<V*> v) {
        // ASSUME v.size()==1
        V *vr = new V;
        vr->order = v[0]->order-1;

        // this probably doesn't work for non-singletons, need true equality check
        for(auto &[k,v] : v[0]->tab) {
            for(auto &[k1,v1] : k->tab) {
                vr->tab[k1] += v*v1;
            }
        }
    }
};

template<class P, class V> FinStochMor<P,V> *FinStochObj<P,V>::samp() {
    if(type.size()!=1) return NULL; // only know how to sample single RVs
    if(type[0]<1) return NULL;

    FinStochMorSamp<P,V> *ret = new FinStochMorSamp<P,V>;
    ret->src=this;
    ret->tgt=new FinStochObj<P,V>;
    ret->tgt->type=type;
    ret->tgt->type[0]--;

    return ret;
}

/* sharpen */
template<class P, class V> struct FinStochMorSharp : public FinStochMor<P,V> {
    using Mor = FinStochMor<P,V>;
    Mor *minner;
    virtual std::vector<V*> eval_on(std::vector<V*> v) {
        std::vector<V*> d = minner->eval_on(v);
        V *vr = new V;
        vr->order = v[0]->order+1;

        vr->tab[v[0]] = 1.0;

        return std::vector<V*> { vr };
    }
};

template<class P, class V> FinStochMor<P,V> *FinStochMor<P,V>::sharp() {
    if(tgt->type.size()!=1) return NULL; // only know how to sample single RVs

    FinStochMorSharp<P,V> *ret = new FinStochMorSharp<P,V>;
    ret->src = src;
    ret->tgt=new FinStochObj<P,V>;
    ret->tgt->type = tgt->type;
    ret->tgt->type[0]++;
    ret->minner = this;
    return ret; 
}

/* condition */
template<class P, class V> struct FinStochMorCond : public FinStochMor<P,V> {
    using Mor = FinStochMor<P,V>;
    Mor *minner;
    std::set<int> os;
    virtual std::vector<V*> eval_on(std::vector<V*> v) {
        V *vr = new V;

        return std::vector<V*> { vr };
    }
};

template<class P, class V> FinStochMor<P,V> *FinStochMor<P,V>::cond(std::set<int> os) {
    FinStochMorCond<P,V> *ret = new FinStochMorCond<P,V>;
    ret->src=new FinStochObj<P,V>;
    ret->src->type = src->type;
    ret->tgt=new FinStochObj<P,V>;
    ret->tgt->type = tgt->type;
    // meanings of indices shift as we erase
    int shift=0;
    for(int i : os) {
        // make this output an input
        ret->src->type.push_back(ret->tgt->type[i-shift]);
        // erase it from the list of outputs
        ret->tgt->type.erase(ret->tgt->type.begin()+(i-shift));
        ++shift;
    }
    ret->minner = this;
    return ret; 
}



