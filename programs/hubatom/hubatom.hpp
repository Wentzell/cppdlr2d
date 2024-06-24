/*!
 *\file hubatom.hpp
 *\brief Header file for Hubbard atom example
 */

#include <cppdlr/cppdlr.hpp>

/*!
 * \brief Driver function for Hubbard atom example, all functions
 *
 * This function expands the Green's function, three-point correlators, and
 * vertex functions in all channels for the Hubbard atom in the DLR. It then
 * computes the polarization function in all channels. It measures the error of
 * all of these representations.
 *
 * \param[in] beta          Inverse temperature
 * \param[in] u             Hubbard interaction
 * \param[in] lambda        DLR cutoff
 * \param[in] eps           Error tolerance
 * \param[in] niomtst       # Matsubara freqs per dim in test grid
 * \param[in] nbos_tst      # bosonic Matsubara freqs in test grid for
 * polarization
 * \param[in] reduced       Use reduced fine Matsubara freq grid
 * \param[in] compressbasis Recompress 2D DLR basis
 * \param[in] niom_dense    # Matsubara freqs in fine grid (only used if
 * reduced=false)
 *
 * \return Vector containing problem parameters and errors, for analysis and
 * plotting
 */
nda::vector<double> hubatom_allfuncs(double beta, double u, double lambda,
                                     double eps, int niomtst, int nbos_tst,
                                     bool reduced, bool compressbasis,
                                     int niom_dense = 0);

/*!
 * \brief Driver function for Hubbard atom example, singlet vertex only
 *
 * This function expands the singlet vertex of the Hubbard atom in the
 * DLR, computes the singlet polarization, and measures their error.
 *
 * \param[in] beta          Inverse temperature
 * \param[in] u             Hubbard interaction
 * \param[in] lambda        DLR cutoff
 * \param[in] eps           Error tolerance
 * \param[in] niomtst       # Matsubara freqs per dim in test grid
 * \param[in] nbos_tst      # bosonic Matsubara freqs in test grid for
 * polarization
 */

void hubatom_lam_s(double beta, double u, double lambda, double eps,
                  int niomtst, int nbos_tst);

// nda::vector<double> hubatom_lam_s(double beta, double u, double lambda,
//                                    double eps, int niomtst, int nbos_tst);

void hubatom_chi_d(double beta, double u, double lambda,
                                     double eps, int niomtst, int nbos_tst);
/*!
 * \defgroup HubSolns
 * \brief Analytical solutions for the Hubbard atom example
 * \note For now, unfortunately, nu means i*nu; change this in the future
 * @{
 */

/*!
 * \brief Green's function
 */
std::complex<double> g_fun(double u, std::complex<double> nu);

/*!
 * \brief Three-point correlator, singlet channel
 */
std::complex<double> chi_s_fun(double u, double beta, std::complex<double> nu1,
                               std::complex<double> nu2);

/*!
 * \brief Three-point correlator, particle-hole channel
 *
 * \param channel +1 (density) or -1 (magnetic)
 */
std::complex<double> chi_ph_fun(double u, double beta, std::complex<double> nu1,
                                std::complex<double> nu2, int channel);

/*!
 * \brief Three-point correlator, density channel
 */
std::complex<double> chi_d_fun(double u, double beta, std::complex<double> nu1,
                               std::complex<double> nu2);

/*!
 * \brief Three-point correlator, magnetic channel
 */
std::complex<double> chi_m_fun(double u, double beta, std::complex<double> nu1,
                               std::complex<double> nu2);

/*!
 * \brief Vertex function, singlet channel
 */
std::complex<double> lam_s_fun(double u, double beta, std::complex<double> nu1,
                               std::complex<double> nu2);

/*!
 * \brief Vertex function, particle-hole channel
 *
 * \param channel +1 (density) or -1 (magnetic)
 */
std::complex<double> lam_ph_fun(double u, double beta, std::complex<double> nu1,
                                std::complex<double> nu2, int channel);

/*!
 * \brief Vertex function, density channel
 */
std::complex<double> lam_d_fun(double u, double beta, std::complex<double> nu1,
                               std::complex<double> nu2);

/*!
 * \brief Vertex function, magnetic channel
 */
std::complex<double> lam_m_fun(double u, double beta, std::complex<double> nu1,
                               std::complex<double> nu2);

/** @} */ // end of HubSolns group
