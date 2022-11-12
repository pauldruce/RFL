#ifndef CLIFFORD_HPP
#define CLIFFORD_HPP

#include <armadillo>
#include <vector>

class Cliff {
public:


    // ============== CONSTRUCTORS, ASSIGNMENT, DESTRUCTOR
    explicit Cliff(int mode);

    Cliff(int p, int q);

    Cliff(const Cliff &C);

    Cliff &operator=(const Cliff &C);

    ~Cliff() = default;;
    // ============== CONSTRUCTORS, ASSIGNMENT, DESTRUCTOR


    // ============== OPERATORS
    Cliff &operator*=(const Cliff &C);

    friend Cliff operator*(Cliff C1, const Cliff &C2) {
        C1 *= C2;
        return C1;
    }
    // ============== OPERATORS


    // ============== GET METHODS
    [[nodiscard]] int get_p() const { return p; }

    [[nodiscard]] int get_q() const { return q; }

    [[nodiscard]] int get_dim_gamma() const { return dim_gamma; }

    [[nodiscard]] arma::cx_mat get_gamma(int i) const { return gamma.at(i); }

    [[nodiscard]] std::vector<arma::cx_mat> get_gamma() const { return gamma; }

    [[nodiscard]] arma::cx_mat get_chiral() const { return chiral; }
    // ============== GET METHODS


    // ============== OTHER METHODS
    void sort_gamma();
    // ============== OTHER METHODS



private:

    int p;
    int q;

    int dim_gamma;

    std::vector<arma::cx_mat> gamma;
    arma::cx_mat chiral;

    void init_gamma();

};


std::ostream &operator<<(std::ostream &out, const Cliff &C);

//void decomp(int p, int q, int &dec);

bool hermiticity(const arma::cx_mat &M1, const arma::cx_mat &M2);


#endif

