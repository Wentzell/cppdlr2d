#include "../src/dlr2d_imfreq.hpp"

#include <fmt/format.h>

using namespace dlr2d;

//TEST(dlr_imfreq, h5_rw) {
int main() {

  double lambda  = 10;      // DLR cutoff
  double eps     = 1e-5;    // DLR tolerance
  auto statistic = Fermion; // Fermionic Green's function

  // Get DLR imaginary frequency object
  auto ifops = imfreq_ops_2d(lambda, eps);

  auto filename = "data_imfreq_ops_h5_rw.h5";
  auto name     = "ifops";

  PRINT(ifops.rank());

  {
    h5::file file(filename, 'w');
    h5::write(file, name, ifops);
  }

  imfreq_ops_2d ifops_ref;
  {
    h5::file file(filename, 'r');
    h5::read(file, name, ifops_ref);
  }

  //// Check equal
  //EXPECT_EQ(ifops.lambda(), ifops_ref.lambda());
  //EXPECT_EQ(ifops.rank(), ifops_ref.rank());
  //EXPECT_EQ_ARRAY(ifops.get_rfnodes(), ifops_ref.get_rfnodes());
  //EXPECT_EQ_ARRAY(ifops.get_ifnodes(), ifops_ref.get_ifnodes());
  //EXPECT_EQ_ARRAY(ifops.get_cf2if(), ifops_ref.get_cf2if());
  //EXPECT_EQ_ARRAY(ifops.get_if2cf_lu(), ifops_ref.get_if2cf_lu());
  //EXPECT_EQ_ARRAY(ifops.get_if2cf_piv(), ifops_ref.get_if2cf_piv());
}
