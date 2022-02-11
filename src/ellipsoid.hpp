/**
   
   Modified version of code by Bojan Nikolic <bojan@bnikolic.co.uk>, Initial version 2010
   
   It has been modified to be used in a stand-alone tool for bounding elipsoid computation
   by Ruediger Ehlers, firstname.lastname@tu-clausthal.de

   The original version and this version are licensed under GNU General Public License version 2

   Computation and use of ellipsoids releated to sets of
   points. References are to Todd and Yildirim, "On Khachiyan's
   Algorithm for the Computation of Minimum Volume Enclosing
   Ellipsoids", 2005

*/
#ifndef _ELIIPSOID_HPP__
#define _ELIIPSOID_HPP__

#include <set>
#undef BOOST_UBLAS_TYPE_CHECK
#define BOOST_UBLAS_TYPE_CHECK 0
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace ublas=boost::numeric::ublas;

/**  Lift to a higher dimensional space that so that points are
   centrally symmetric.

   The input data are arranged in columns of A

   c.f. Eq 9
*/
void Lift(const boost::numeric::ublas::matrix<double> &A,
        boost::numeric::ublas::matrix<double> &Ap);


/** Compute the result of operator Lambda

  c.f. eq. 13

  \param Lambdap The result is stored here
*/
void KaLambda(const boost::numeric::ublas::matrix<double> &Ap,
            const boost::numeric::ublas::vector<double> &p,
            ublas::matrix<double> &Lambdap);

/** Calculate the inverse of a result of lambda operator
*/
void InvertLP(const boost::numeric::ublas::matrix<double> &Lambdap,
            boost::numeric::ublas::matrix<double> &LpInv);


/** Iterate the Khachiyan algorithm one step

  \param Ap is the lifted data set

  \param p is the current guess, and the iterated guess is stored
  here

  \returns Error estimate
*/
double KhachiyanIter(const boost::numeric::ublas::matrix<double> &Ap,
                   boost::numeric::ublas::vector<double> &p);


/** Invert from the dual problem back to real space

  \param Q The ellipsiod matrix
  \param c The centre of the ellipse
*/
void KaInvertDual(const boost::numeric::ublas::matrix<double> &A,
                const boost::numeric::ublas::vector<double> &p,
                boost::numeric::ublas::matrix<double> &Q,
                boost::numeric::ublas::vector<double> &c
                );


/** \param Compute minimum containing ellipsoid using the Khachiyan
  algorithm

  \param A Input point set as columns of the matrix

  \param eps Tolerance of the estimated ellipsoid: it's volume
  will be approx (1+eps) of the true smallest ellipsoid

  \param maxiter Maximum number of iterations to make before
  stopping the algorithm

  \param Q The matrix describing the axes of the resulting
  ellipsoid

  \param c The centre of the resulting ellipsoid
*/
double KhachiyanAlgo(const ublas::matrix<double> &A,
                   double eps,
                   size_t maxiter,
                   boost::numeric::ublas::matrix<double> &Q,
                   boost::numeric::ublas::vector<double> &c);


/** \brief Convenience class to store the results of minimum
  covering ellipsoid calculation
*/
struct KhachiyanEllipsoid
{
    boost::numeric::ublas::matrix<double> Q;
    boost::numeric::ublas::vector<double> c;
};



#endif
