#include "tab2D.hxx"
#include "many_tab2D.hxx"

template class Table2D<ELEM_TYPE>;

template class Table2D<BIG_ELEM_TYPE>;
template class CManyTab2D<BIG_ELEM_TYPE>;
template class Table2D<double>;
template class Table2D<short>;
template class Table2D<float>;
template class Table2D<unsigned char>;
