#include<vector>
#include<map>
#include<unordered_map>
#include<functional>
#include<set>
#include<armadillo>

/* abstract interface to markov categories */
struct Mor;
struct Obj {
    virtual Mor* id()=0;
    virtual int tensor_length()=0;
    virtual Mor* copy_n(int n)=0;
    virtual bool is_equal(Obj *other)=0;
    virtual Obj* tensor(Obj *other)=0;
    virtual Mor* term()=0; // morphism to terminal object

    virtual Mor* select(int n)=0; // morphism to nth tensee, technically derivable but we chose to not expose the tensor structure at this stage...

    virtual void debug_print()=0;

    // derived from the others
    Mor* permute_project(std::vector<int> v);
};
struct Mor {
    virtual Obj* dom()=0;
    virtual Obj* cod()=0;
    virtual Mor* cond(std::vector<Obj*> os)=0;
    virtual Mor* compose(Mor* other)=0;
    virtual Mor* tensor(Mor *other)=0;

    virtual void debug_print()=0;
};


Mor *Obj::permute_project(std::vector<int> v)
{
    Mor *delta = copy_n(v.size());
    Mor *join = NULL;
    for(int i=0;i<v.size();++i) {
        Mor *selectone = select(v[i]);
        
        if(join) join = join->tensor(selectone);
        else join = selectone;
    }
    return delta->compose(join);
}


/* an enumerable domain for finstoch */
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
        printf("FinStochObj %X: ", this);
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
    virtual FinStochMor* cond(std::vector<Obj*> os) {};
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
        printf("== FinStochMor %X ==\n", this);
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
