#include "mcat.h"

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

    return 0;
}
