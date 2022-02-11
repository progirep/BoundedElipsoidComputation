/**
   Modified version of code by Bojan Nikolic <bojan@bnikolic.co.uk>, Initial version 2010
   
   It has been modified to be used in a stand-alone tool for bounding elipsoid computation
   by Ruediger Ehlers, firstname.lastname@tu-clausthal.de

   The original version and this version are licensed under GNU General Public License version 2
   
   Computation and use of ellipsoids releated to sets of points
*/

#undef BOOST_UBLAS_TYPE_CHECK
#define BOOST_UBLAS_TYPE_CHECK 0
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <iostream>
#include <boost/numeric/ublas/io.hpp>

#include "ellipsoid.hpp"
//#include "../mcpoint.hpp"


template<class T>
bool InvertMatrix(const ublas::matrix<T> &input,
                ublas::matrix<T> &inverse)
{
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;
    matrix<T> A(input);
    pmatrix pm(A.size1());
    int res = lu_factorize(A,pm);
    if( res != 0 ) return false;
    inverse.assign(ublas::identity_matrix<T>(A.size1()));
    lu_substitute(A, pm, inverse);
    return true;
}

void InvertLP(const ublas::matrix<double> &Lambdap,
            ublas::matrix<double> &LpInv)
{
    bool res=InvertMatrix(Lambdap, LpInv);
    if (not res)
    {
     // throw an error of your choice here!
     // throw MatrixErr("Could not invert matrix: ",
     //                 Lambdap);
    }
}

void Lift(const ublas::matrix<double> &A,
        ublas::matrix<double> &Ap)
{
    Ap.resize(A.size1()+1,
              A.size2());
    ublas::matrix_range<ublas::matrix<double> >
      sub(Ap,
          ublas::range(0, A.size1()),
          ublas::range(0, A.size2()));
    sub.assign(A);
    ublas::row(Ap, Ap.size1()-1)=ublas::scalar_vector<double>(A.size2(),1.0);

}

void genDiag(const ublas::vector<double> &p,
           ublas::matrix<double> &res)
{
    res.assign(ublas::zero_matrix<double>(p.size(),
                                          p.size()));
    for(size_t i=0; i<p.size(); ++i)
    {
      res(i,i)=p(i);
    }
    }

void KaLambda(const ublas::matrix<double> &Ap,
            const ublas::vector<double> &p,
            ublas::matrix<double> &Lambdap)
{

    ublas::matrix<double> dp(p.size(), p.size());
    genDiag(p, dp);

    dp=ublas::prod(dp, ublas::trans(Ap));
    Lambdap=ublas::prod(Ap,
                        dp);
}

double KhachiyanIter(const ublas::matrix<double> &Ap,
                   ublas::vector<double> &p)
{
    /// Dimensionality of the problem
    const size_t d=Ap.size1()-1;

    ublas::matrix<double> Lp;
    ublas::matrix<double> M;
    KaLambda(Ap, p, Lp);
    ublas::matrix<double> ILp(Lp.size1(), Lp.size2());
    InvertLP(Lp, ILp);
    M=ublas::prod(ILp, Ap);
    M=ublas::prod(ublas::trans(Ap), M);

    double maxval=0;
    size_t maxi=0;
    for(size_t i=0; i<M.size1(); ++i)
    {
      if (M(i,i) > maxval)
      {
        maxval=M(i,i);
        maxi=i;
      }
    }
    const double step_size=(maxval -d - 1)/((d+1)*(maxval-1));
    ublas::vector<double> newp=p*(1-step_size);
    newp(maxi) += step_size;

    const double err= ublas::norm_2(newp-p);
    p=newp;
    return err;

}

void KaInvertDual(const ublas::matrix<double> &A,
                const ublas::vector<double> &p,
                ublas::matrix<double> &Q,
                ublas::vector<double> &c
                )
{
    const size_t d=A.size1();
    ublas::matrix<double> dp(p.size(), p.size());
    genDiag(p, dp);

    ublas::matrix<double> PN=ublas::prod(dp, ublas::trans(A));
    PN=ublas::prod(A, PN);

    ublas::vector<double> M2=ublas::prod(A, p);
    ublas::matrix<double> M3=ublas::outer_prod(M2, M2);

    ublas::matrix<double> invert(PN.size1(), PN.size2());
    InvertLP(PN- M3, invert);

    Q.assign( 1.0/d *invert);
    c=ublas::prod(A, p);


}




