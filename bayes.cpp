#include "mcat.h"
#include "finstoch.h"

Mor* ibf(Mor *f0, Mor* f, Mor *g, std::vector<Mor*> ys)
{
    if(ys.size()==1) {
        // base case
        Mor *y1 = ys.back();
        Mor *inner =  f0->compose(
                        f0->cod()->copy_n(2)->compose(
                           f0->cod()->id()->tensor(g) ));
        Mor *cond = inner->cond(1);
        return y1->compose(cond);
    } else {
        Mor *yn = ys.back();
        ys.pop_back();
        Mor *prev = ibf(f0,f,g,ys);
        prev->debug_print();
        Mor *inner = prev->compose(
                       f->compose(
                          f->cod()->copy_n(2)->compose(
                             f->cod()->id()->tensor(g) )));
        Mor *cond = inner->cond(1);
        return yn->compose(cond);
    }
}

int main(int argc, char* argv[])
{
    FinStochObj *two = new FinStochObj(2);
    FinStochMor *f = new FinStochMor();
    f->o_dom = two;
    f->o_cod = two;
    f->data = { {0.8,0.2}, {0.0,1.0} };
    f->debug_print();

    FinStochMor *g = new FinStochMor();
    g->o_dom = two;
    g->o_cod = two;
    g->data = { {0.9,0.1}, {0.1,0.9} };
    g->debug_print();

    FinStochMor *f0 = two->det_const({0});
    f0->debug_print();

    srand((unsigned)time(NULL));
    std::vector<Mor*> states;
    int st=0,stf;
    for(int i=0;i<10;++i) {
        if(st==0) {
            if(!(rand()%5)) st=1;
        }
        if(rand()%10) stf=st; else stf=!st;
        states.push_back(two->det_const({stf}));
        printf("%d",stf);
    }
    printf("\n");

    FinStochMor *m = (FinStochMor*)ibf(f0,f,g,states);
    //m->debug_print();

    return 0;
}
