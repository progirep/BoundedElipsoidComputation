#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <signal.h>
#include <fstream>
#include <limits>
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


void printVolume(ublas::matrix<double> Q /* Copied on purpose as it's chanced in-place*/, const ublas::vector<double> &c, const ublas::matrix<double> *A) {

    // Correct matrix
    double maxDistance = 0.0;
    unsigned int nofDimensions = Q.size1();
    for (unsigned int i=0;i<A->size2();i++) {
        // Compute
        double result = 0.0;
        for (unsigned int j=0;j<nofDimensions;j++) {
            for (unsigned int m=0;m<nofDimensions;m++) {
                result += Q(j,m)*( (*A)(j,i) - c(j))*( (*A)(m,i) - c(m));
            }
        }
        maxDistance = std::max(maxDistance,result);
    }

    // Scale matrix
    for (unsigned int y=0;y<nofDimensions;y++) {
        for (unsigned int x=0;x<nofDimensions;x++) {
            Q(x,y) /= maxDistance;
        }
    }

    ublas::permutation_matrix<size_t> permutationMatrix(Q.size1());
    if(ublas::lu_factorize(Q,permutationMatrix)) {
        std::cerr << "0\n";
    } else {
        // Source for sphere volume: https://math.stackexchange.com/questions/332391/volume-of-hyperellipsoid/332434
        double nD2 = Q.size1()/2.0;
        double volume = nD2*pow(M_PI,nD2)/tgamma(nD2);
        for(unsigned int i=0; i<Q.size1(); i++) {
            //std::cout << "[" << sqrt(abs(Q(i,i))) << "]";
            volume *= 1.0/sqrt(fabs(Q(i,i)));
        }
        std::cout << std::abs(volume) << "\n";
    }
}

/* Signal handler */
volatile sig_atomic_t notAborted = true;
void sigintHandler(int) {
    notAborted = false;
}

double KhachiyanAlgo(const ublas::matrix<double> &A,
                   double eps,
                   size_t maxiter,
                   ublas::matrix<double> &Q,
                   ublas::vector<double> &c,
                   int printVolumeEvery)
{
    ublas::vector<double> p=ublas::scalar_vector<double>(A.size2(), 1.0)*(1.0/A.size2());

    ublas::matrix<double> Ap;
    Lift(A, Ap);

    double ceps=std::numeric_limits<double>::max();
    for (size_t i=1;  i<=maxiter && (ceps>eps) && notAborted; ++i) {
      ceps=KhachiyanIter(Ap, p);
      if ((printVolumeEvery>-1) && (i % printVolumeEvery)==0) {
        std::cerr << "Current Volume: ";
        KaInvertDual(A, p, Q, c);
        printVolume(Q,c,&A);
      }
    }
    if (!notAborted) {
        std::cerr << "Note: Aborted computation due to SIGINT.\n";
    }

    KaInvertDual(A, p, Q, c);

    return ceps;


}


void readMatrixFromFile(ublas::matrix<double> **A, unsigned int *nofDimensions, std::istream *is, uint32_t nofInitialLinesToBeSkipped) {


    // Read from Stdin - Unknown length of input. Need to buffer it.
    std::list<std::vector<double> > bufferedLines;
    int nofLine = 1;
    std::string nextLine;

    for (uint32_t i=0;i<nofInitialLinesToBeSkipped;i++) {
        std::string wasted;
        if (!std::getline(*is,wasted)) {
            std::cerr << "Error: The input does not have as many lines as they are to be skipped.\n";
            std::exit(1);
        }
    }

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

    double epsilon = 0.00001;
    int maxNofIterations = 1000;
    int printVolumeEveryNthIteration = -1;
    bool installSIGINTHandler = true;
    uint32_t nofInitialLinesToBeSkipped = 0;

    std::string inFileName = ""; // empty string represents: No file name given.
    for (int i=1;i<nofArgs;i++) {
        std::string thisArg = args[i];
        if (thisArg.size()>0) {
            if (thisArg[0]=='-') { // Ignore empty parameters

                // Option Parameter

                if (thisArg=="--maxPrecision") {
                    if (i>=nofArgs-1) {
                        std::cerr << "Error: Expecting parameter after '--maxPrecision'\n";
                        return 1;
                    }
                    std::istringstream is(args[i+1]);
                    is >> epsilon;
                    if (is.fail()) {
                        std::cerr << "Error: Could not parse floating point number '" << args[i+1] << "'\n";
                        return 1;
                    }
                    i++;
                    if (epsilon<0) {
                        std::cerr << "Note: Running a fixed number of iterations due to negative epsilon value.\n";
                    }
                } else if (thisArg=="--maxIterations") {
                    if (i>=nofArgs-1) {
                        std::cerr << "Error: Expecting parameter after '--maxIterations'\n";
                        return 1;
                    }
                    std::istringstream is(args[i+1]);
                    is >> maxNofIterations;
                    if (is.fail()) {
                        std::cerr << "Error: Could not parse integer number '" << args[i+1] << "'\n";
                        return 1;
                    }
                    i++;
                    if (maxNofIterations<1) {
                        std::cerr << "Error: Maximum number of iterations needs to be positive.\n";
                        return 1;
                    }

                } else if (thisArg=="--skipfirstNLines") {
                    if (i>=nofArgs-1) {
                        std::cerr << "Error: Expecting parameter after '--skipfirstNLines'\n";
                        return 1;
                    }
                    std::istringstream is(args[i+1]);
                    is >> nofInitialLinesToBeSkipped;
                    if (is.fail()) {
                        std::cerr << "Error: Could not parse integer number '" << args[i+1] << "'\n";
                        return 1;
                    }
                    i++;

                } else if (thisArg=="--printVolumeEvery") {
                    if (i>=nofArgs-1) {
                        std::cerr << "Error: Expecting parameter after '--printVolumeEvery'\n";
                        return 1;
                    }
                    std::istringstream is(args[i+1]);
                    is >> printVolumeEveryNthIteration;
                    if (is.fail()) {
                        std::cerr << "Error: Could not parse integer number '" << args[i+1] << "'\n";
                        return 1;
                    }
                    i++;
                    if (printVolumeEveryNthIteration<1) {
                        std::cerr << "Error: The number n for which the elipsoid volume is printed every n iterations needs to be positive.\n";
                        return 1;
                    }

                } else if (thisArg=="--doNotInstallSIGINTHandler") {
                    installSIGINTHandler = false;
                } else {

                    // Option not recognized.
                    std::cerr << "Error: Did not understand parameter: '" << thisArg << "'.\n";
                    return 1;
                }

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

    // SIGINT Handler
    struct sigaction act;
    if (installSIGINTHandler) {
        act.sa_handler = sigintHandler;
        sigaction(SIGINT, &act, NULL);
    }

    // Input Matrix
    ublas::matrix<double> *A;
    unsigned int nofDimensions = 0;

    if (inFileName=="") {
        readMatrixFromFile(&A,&nofDimensions,&(std::cin),nofInitialLinesToBeSkipped);
    } else {
        std::ifstream inFile(inFileName);
        readMatrixFromFile(&A,&nofDimensions,&(inFile),nofInitialLinesToBeSkipped);
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
                                     Q, c,printVolumeEveryNthIteration);

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

    std::cout << "Result - Volume Elipsoid: ";
    ublas::permutation_matrix<size_t> permutationMatrix(Q.size1());
    if(ublas::lu_factorize(Q,permutationMatrix)) {
        std::cout << "0\n";
    } else {
        // Source for sphere volume: https://math.stackexchange.com/questions/332391/volume-of-hyperellipsoid/332434
        double nD2 = Q.size1()/2.0;
        double volume = nD2*pow(M_PI,nD2)/tgamma(nD2);
        for(unsigned int i=0; i<Q.size1(); i++) {
            volume *= 1.0/sqrt(fabs(Q(i,i)));
        }
        std::cout << std::abs(volume) << "\n";
    }


    delete A;

}
