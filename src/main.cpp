#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <fstream>
#include "ellipsoid.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
/*
 * Tool for computing an elipsoid containing a set of n-dimensional points
 * from a given point list.
 *
 * The points are either read from a file or via stdin.
 *
 * Licensed under GPL v2. See "LICENSE.txt" for details
 */


double KhachiyanAlgo(const ublas::matrix<double> &A,
                   double eps,
                   size_t maxiter,
                   ublas::matrix<double> &Q,
                   ublas::vector<double> &c)
{
    ublas::vector<double> p=ublas::scalar_vector<double>(A.size2(), 1.0)*(1.0/A.size2());

    ublas::matrix<double> Ap;
    Lift(A, Ap);

    double ceps=eps*2;
    for (size_t i=0;  i<maxiter && ceps>eps; ++i) {
      ceps=KhachiyanIter(Ap, p);
    }

    KaInvertDual(A, p, Q, c);

    return ceps;


}


void readMatrixFromFile(ublas::matrix<double> **A, unsigned int *nofDimensions, std::istream *is) {


    // Read from Stdin - Unknown length of input. Need to buffer it.
    std::list<std::vector<double> > bufferedLines;
    int nofLine = 1;
    std::string nextLine;
    while (std::getline(*is,nextLine)) {
        // Interpret line.

        std::istringstream lineReader(nextLine);
        std::vector<double> thisOne;
        while (!(lineReader.eof())) {
            double a;
            while (lineReader >> a) {
                thisOne.push_back(a);
            }
            if (lineReader.bad()) {
                std::cerr << "Error parsing input line no. " << nofLine << ": '" << nextLine << "'\n";
                std::exit(1);
            }
        }
        if (thisOne.size()>0) {
            if (*nofDimensions!=0) {
                if (thisOne.size()!=*nofDimensions) {
                    std::cerr << "Error: Number of dimensions in line no. " << nofLine << " is not consistent to the first line of the data.\n";
                    std::exit(1);
                }
            } else {
                *nofDimensions = thisOne.size();
            }
            bufferedLines.push_back(thisOne);
        }
        nofLine++;
    }

    // Compute A
    *A = new ublas::matrix<double>(*nofDimensions,bufferedLines.size());

    size_t j=0;
    for (std::list<std::vector<double> >::const_iterator i=bufferedLines.begin();
         i != bufferedLines.end();++i) {
      for(size_t k=0; k <*nofDimensions ;++k)
        (**A)(k,j)=(*i)[k];
      ++j;
    }

}




int main(int nofArgs, char **args) {

    double epsilon = 0.00000001;
    int maxNofIterations = 1000000;

    std::string inFileName = ""; // empty string represents: No file name given.
    for (int i=1;i<nofArgs;i++) {
        std::string thisArg = args[i];
        if (thisArg.size()>0) {
            if (thisArg[0]=='-') { // Ignore empty parameters

                // Option Parameter
                std::cerr << "Error: Did not understand parameter: '" << thisArg << "'.\n";
                return 1;

            } else {
                // File name
                if (inFileName=="") {
                    inFileName = thisArg;
                } else {
                    std::cerr << "Error: Multiple input file names given.\n";
                    return 1;
                }
            }
        }
    }

    // Input Matrix
    ublas::matrix<double> *A;
    unsigned int nofDimensions = 0;

    if (inFileName=="") {
        readMatrixFromFile(&A,&nofDimensions,&(std::cin));
    } else {
        std::ifstream inFile(inFileName);
        readMatrixFromFile(&A,&nofDimensions,&(inFile));
    }


    ublas::matrix<double> Q(nofDimensions,nofDimensions);
    ublas::vector<double> c(nofDimensions);

    std::cout << "Input - Matrix:\n";
    for (unsigned int y=0;y<A->size2();y++) {
        for (unsigned int x=0;x<nofDimensions;x++) {
            if (x!=0) std::cout << " ";
            std::cout << (*A)(x,y);
        }
        std::cout << "\n";
    }

    const double ceps=KhachiyanAlgo(*A, epsilon, maxNofIterations,
                                     Q, c);

    // Result is in Q and c, ceps is the actual epsilon
    // Print matrix
    std::cout << "Result - Center vector:\n";
    for (unsigned int x=0;x<nofDimensions;x++) {
        if (x!=0) std::cout << " ";
        std::cout << c(x);
    }
    std::cout << "\n";

    std::cout << "Testing Elipsoid:\n";
    double maxDistance = 0.0;
    for (unsigned int i=0;i<A->size2();i++) {
        for (unsigned int x=0;x<nofDimensions;x++) {
            std::cout << (*A)(x,i) << " ";
        }
        std::cout << ": ";

        // Compute
        double result = 0.0;
        for (unsigned int j=0;j<nofDimensions;j++) {
            for (unsigned int m=0;m<nofDimensions;m++) {
                result += Q(j,m)*( (*A)(j,i) - c(j))*( (*A)(m,i) - c(m));
            }
        }
        maxDistance = std::max(maxDistance,result);

        std::cout << result << "\n";
    }

    // Scale matrix
    for (unsigned int y=0;y<nofDimensions;y++) {
        for (unsigned int x=0;x<nofDimensions;x++) {
            Q(x,y) /= maxDistance;
        }
    }

    std::cout << "Result - Matrix:\n";
    for (unsigned int y=0;y<nofDimensions;y++) {
        for (unsigned int x=0;x<nofDimensions;x++) {
            if (x!=0) std::cout << " ";
            std::cout << Q(y,x);
        }
        std::cout << "\n";
    }

    std::cout << "Result - CEps: " << ceps << "\n";

    std::cout << "Result - Volume Elipsoid:\n";
    ublas::permutation_matrix<size_t> permutationMatrix(Q.size1());
    if(ublas::lu_factorize(Q,permutationMatrix)) {
        std::cout << "Error: Elipsoid has no volume.\n";
        return 1;
    } else {
        // Source for sphere volume: https://math.stackexchange.com/questions/332391/volume-of-hyperellipsoid/332434
        double nD2 = nofDimensions/2.0;
        double volume = nD2*pow(M_PI,nD2)/tgamma(nD2);
        for(unsigned int i=0; i<Q.size1(); i++) {
            volume *= 1.0/Q(i,i);
        }
        std::cout << std::abs(volume) << "\n";
    }

    delete A;

}
