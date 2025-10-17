#include "ElGamal.hpp"
#include "cramer_shoup.hpp"
#include "hashed_ElGamal.hpp"
#include "twin_hashed_ElGamal.hpp"

int main() {
  test_ElGamal();
  test_cramer_shoup();
  test_hashed_ElGamal();
  test_twin_hashed_ElGamal();
  return 0;
}