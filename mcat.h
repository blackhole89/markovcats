#ifndef _MCAT_H_
#define _MCAT_H_

#include<vector>

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
    virtual Mor* cond(int n)=0; //switch n last tensees of cod into dom
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

#endif


