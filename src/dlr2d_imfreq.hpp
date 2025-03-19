#pragma once

#include <nda/nda.hpp>

#include "dlr2d.hpp"

namespace dlr2d {

  /**
  * @class imfreq_ops
  * @brief Class responsible for all 2D DLR imaginary frequency operations,
  * including building imaginary frequency grid and transformations.
  *
  * \note First dimension of all Green's function and coefficient arrays must be
  * DLR rank r.
  */

  class imfreq_ops_2d {

    public:
    /** 
    * @brief Constructor for imfreq_ops
    * 
    * @param[in] lambda DLR cutoff parameter
    * @param[in] eps Error tolerance
    */
    imfreq_ops_2d(double lambda, double eps) : lambda_(lambda), eps_(eps) {

      dlr2d_if = build_dlr2d_if(lambda, eps);

      r      = dlr2d_if.size();
      dlr_rf = cppdlr::build_dlr_rf(lambda, eps, false);
      cf2if  = build_cf2if(lambda, dlr_rf, dlr2d_if);

      if2cf.lu  = fmatrix(cf2if);
      if2cf.piv = nda::vector<int>(r);
      nda::lapack::getrf(if2cf.lu, if2cf.piv);
    }

    imfreq_ops_2d(double lambda, double eps, nda::vector_const_view<double> dlr_rf, //
                  nda::array<int, 2> dlr2d_if,                                      //
                  fmatrix_const_view cf2if,                                         //
                  fmatrix_const_view if2cf_lu,                                      //
                  nda::vector_const_view<int> if2cf_piv)
       : lambda_(lambda), r(cf2if.extent(1)), dlr_rf(dlr_rf), dlr2d_if(dlr2d_if), cf2if(cf2if), if2cf{if2cf_lu, if2cf_piv} {};

    imfreq_ops_2d() = default;

    /** 
    * @brief Transform values of Green's function G on 2D DLR imaginary frequency grid to
    * 2D DLR coefficients
    *
    * @param[in] g Values of G on 2D DLR imaginary frequency grid
    *
    * @return 2D DLR coefficients of G
    */
    template <nda::MemoryArray T> typename T::regular_type vals2coefs(T const &g) const {

      // MemoryArray type can be nda vector, matrix, array, or view of any of
      // these; taking a regular_type converts, for example, a matrix view to a
      // matrix.

      if (r != g.shape(0)) throw std::runtime_error("First dim of g != # DLR imaginary frequency nodes.");

      // Make a copy of the data in Fortran Layout as required by getrs
      auto gf = nda::array<nda::get_value_t<T>, nda::get_rank<T>, nda::F_layout>(g);

      // Reshape as matrix_view with r rows
      auto gfv = nda::reshape(gf, r, g.size() / r);

      // Solve linear system (multiple right hand sides) to convert vals -> coeffs
      nda::lapack::getrs(if2cf.lu, gfv, if2cf.piv);

      return gf(nda::range(r), nda::ellipsis());
    }

    /** 
    * @brief Transform 2D DLR coefficients of Green's function G to values on 2D DLR
    * imaginary frequency grid
    *
    * @param[in] gc 2D DLR coefficients of G
    *
    * @return Values of G on 2D DLR imaginary frequency grid
    */
    template <nda::MemoryArray T, nda::Scalar S = nda::get_value_t<T>> make_cplx_t<T> coefs2vals(T const &gc) const {

      if (r != gc.shape(0)) throw std::runtime_error("First dim of gc != DLR rank r.");

      // Reshape gc to a matrix w/ first dimension r
      auto gc_rs = nda::reshape(gc, r, gc.size() / r);

      // Apply coeffs -> vals matrix
      auto g = cf2if * nda::matrix_const_view<S>(gc_rs);

      // Get output shape
      std::array<long, T::rank> shape_out;
      shape_out[0] = r;
      for (int i = 1; i < T::rank; ++i) { shape_out[i] = gc.shape(i); }

      // Reshape to original dimensions and return
      return nda::reshape(g, shape_out);
    }

    /** 
    * @brief Evaluate 2D DLR expansion of G, given by its DLR coefficients, at imaginary
    * frequency point 
    *
    * @param[in] gc    2D DLR coefficients of G
    * @param[in] n, m  Pair of Matsubara Frequencies for the evaluation
    *
    * @return Value of G at @p Matsubara frequency index pair (n, m)
    */
    template <nda::MemoryArray T, nda::Scalar S = nda::get_value_t<T>> auto coefs2eval(T const &gc, int n, int m) const {

      if (r != gc.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Scalar-valued Green's functions are handled differently than matrix-valued Green's functions

      if constexpr (T::rank == 1) {

        // Evaluate DLR expansion FIXME
        auto g = dcomplex{0};

        return g;
      } else {

        // Reshape gc to matrix w/ first dimension r
        auto gc_rs = nda::reshape(gc, r, gc.size() / r);

        // Get output shape
        std::array<long, T::rank - 1> shape_out;
        for (int i = 0; i < T::rank - 1; ++i) { shape_out[i] = gc.shape(i + 1); }

        // Get vector of evaluation of DLR expansion at a point FIXME
        auto g = nda::matrix<S>{};

        // Evaluate DLR expansion, reshape to original dimensions (with first
        // dimension summed out), and return
        return nda::reshape(g, shape_out);
      }
    }

    /** 
    * @brief Get DLR imaginary frequency nodes
    *
    * @return DLR imaginary frequency nodes
    */
    nda::array_const_view<int, 2> get_ifnodes() const { return dlr2d_if; };
    std::pair<int, int> get_ifnodes(int i) const { return {dlr2d_if(i, 0), dlr2d_if(i, 1)}; };

    /**
    * @brief Get DLR real frequency nodes
    *
    * @return DLR real frequency nodes
    */
    nda::vector_const_view<double> get_rfnodes() const { return dlr_rf; };
    double get_rfnodes(int i) const { return dlr_rf(i); };

    /**
    * @brief Get transformation matrix from DLR coefficients to values at DLR imaginary frequency nodes
    *
    * @return Transformation matrix
    */
    fmatrix_const_view get_cf2if() const { return cf2if; };

    /**
    * @brief Get LU factors of transformation matrix from DLR imaginary frequency values to coefficients
    *
    * @return LU factors
    */
    fmatrix_const_view get_if2cf_lu() const { return if2cf.lu; };

    /**
    * @brief Get LU pivots of transformation matrix from DLR imaginary frequency values to coefficients
    *
    * @return LU pivots
    */
    nda::vector_const_view<int> get_if2cf_piv() const { return if2cf.piv; };

    /** 
    * @brief Get DLR rank
    *
    * @return DLR rank
    */
    int rank() const { return r; }
    double lambda() const { return lambda_; }
    double eps() const { return eps_; }

    private:
    double lambda_;              ///< Energy cutoff divided by temperature
    double eps_;                 ///< Error tolerance
    int r;                       ///< DLR rank
    nda::vector<double> dlr_rf;  ///< 1D DLR frequencies
    nda::array<int, 2> dlr2d_if; ///< 2D DLR imaginary frequency nodes
    fmatrix cf2if;               /// Transformation matrix from DLR coefficients to values at DLR imaginary frequency nodes

    /**
    * @brief Struct for transformation from DLR imaginary frequency values to coefficients
    */
    struct {
      fmatrix lu;           ///< LU factors (LAPACK format) of imaginary frequency vals -> coefs matrix
      nda::vector<int> piv; ///< LU pivots (LAPACK format) of imaginary frequency vals -> coefs matrix
    } if2cf;

    // -------------------- serialization -------------------

    public:
    /**
     * Serialize the object into an archive by serializing all its members.
     * The archive parameter must support the operator& to serialize each member.
     *
     * @param[in] ar Archive to serialize into
     */
    void serialize(auto &ar) const { ar & lambda_ & r & dlr_rf & dlr2d_if & cf2if & if2cf.lu & if2cf.piv; }

    /**
     * Deserialize an object from the archive. This will initialize all members.
     * The archive parameter must support the operator& to deserialize the members
     * in the order they were serialized.
     *
     * @param[in] ar Archive to deserialize from
     */
    void deserialize(auto &ar) { ar & lambda_ & r & dlr_rf & dlr2d_if & cf2if & if2cf.lu & if2cf.piv; }

    // -------------------- hdf5 -------------------

    static std::string hdf5_format() { return "cppdlr::imfreq_ops"; }

    friend void h5_write(h5::group fg, std::string const &subgroup_name, imfreq_ops_2d const &m) {

      h5::group gr = fg.create_group(subgroup_name);
      write_hdf5_format_as_string(gr, "cppdlr::imfreq_ops");

      h5::write(gr, "lambda", m.lambda());
      h5::write(gr, "eps", m.eps());
      h5::write(gr, "rf", m.get_rfnodes());
      h5::write(gr, "if", m.get_ifnodes());
      h5::write(gr, "cf2if", m.get_cf2if());
      h5::write(gr, "if2cf_lu", m.get_if2cf_lu());
      h5::write(gr, "if2cf_piv", m.get_if2cf_piv());
    }

    friend void h5_read(h5::group fg, std::string const &subgroup_name, imfreq_ops_2d &m) {

      h5::group gr = fg.open_group(subgroup_name);
      assert_hdf5_format_as_string(gr, "cppdlr::imfreq_ops", true);

      double lambda, eps;
      nda::vector<double> rf_;
      nda::array<int, 2> if_;
      fmatrix cf2if_;
      fmatrix if2cf_lu;
      nda::vector<int> if2cf_piv;

      h5::read(gr, "lambda", lambda);
      h5::read(gr, "eps", eps);
      h5::read(gr, "rf", rf_);
      h5::read(gr, "if", if_);
      h5::read(gr, "cf2if", cf2if_);
      h5::read(gr, "if2cf_lu", if2cf_lu);
      h5::read(gr, "if2cf_piv", if2cf_piv);

      m = imfreq_ops_2d(lambda, eps, rf_, if_, cf2if_, if2cf_lu, if2cf_piv);
    }
  };

} // namespace dlr2d
