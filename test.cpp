#include "mcat.h"
#include "finstoch.h"

int main(int argc, char* argv[])
{
//    FinStoch<float, FinStochVal<float, 2> > fs;

    FinStochObj *two = new FinStochObj(2);
    FinStochObj *three = new FinStochObj(3);
    FinStochObj *tt = two->tensor(three);
    tt->select(1)->debug_print();
    two->copy_n(3)->debug_print();
    tt->permute_project({1,0})->debug_print();
    tt->permute_project({1})->debug_print();

    tt->det_const({0,2})->compose(tt->permute_project({1}))->debug_print();

    FinStochObj *ttt = tt->tensor(three);
    ttt->det_const({0,1,2})->compose(ttt->permute_project({0,2,1}))->compose(ttt->permute_project({1}))->debug_print();

    tt->permute_project({1,0})->cond(1)->debug_print();

    FinStochObj *twosq = two->tensor(two);

    // unif distr on [2]x[2]
    FinStochMor *m_unif_twosq = new FinStochMor();
    m_unif_twosq->o_dom = new FinStochObj(); // terminal
    m_unif_twosq->o_cod = twosq;
    m_unif_twosq->data = { { 0.25, 0.25, 0.25, 0.25 } };


    FinStochMor *m_or = new FinStochMor();
    m_or->o_dom = twosq;
    m_or->o_cod = two;
    m_or->data = { { 1.0, 0.0}, {0.0, 1.0}, {0.0,1.0}, {0.0,1.0} };
    m_or->debug_print();
    FinStochMor *m_preserve = (FinStochMor*)twosq->copy_n(2)->compose(twosq->id()->tensor(m_or));
    m_preserve = m_unif_twosq->compose(m_preserve);
    m_preserve->debug_print();

    FinStochMor *m_cond = m_preserve->cond(1);
    m_cond->debug_print();

    return 0;
}
