#ifndef HASHED_ELGAMAL_H
#define HASHED_ELGAMAL_H

#include "curve/curve.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <pbc/pbc.h>
#include <pbc/pbc_field.h>
#include <pbc/pbc_pairing.h>
#include <string>

void test_hashed_ElGamal() {
  printf(
      "------------------------- hashed ElGamal -------------------------\n");

  /*
   * sys_gen: choose a cyclic group, a generator g and a hash function H
   */
  pairing_t pairing_;
  pbc_param_t par;
  pbc_param_init_set_str(par, a_param.c_str());
  pairing_init_pbc_param(pairing_, par);

  element_t g;
  element_init_G1(g, pairing_);
  element_random(g);

  std::hash<std::string> hash;

  /*
   * key_gen: sample alpha, sk = alpha, pk = g_1=g^alpha
   */
  element_t alpha, g_1;
  element_init_Zr(alpha, pairing_);
  element_init_G1(g_1, pairing_);
  element_random(alpha);
  element_pow_zn(g_1, g, alpha);

  /*
   * encrypt: ct = (C_1 = g^r, C_2 = H(g_1^r) xor m)
   */
  size_t mes = 12345;
  printf("明文：%ld\n", mes);

  element_t r;
  element_init_Zr(r, pairing_);
  element_random(r);

  // compute C_1
  element_t C_1;
  element_init_G1(C_1, pairing_);
  element_pow_zn(C_1, g, r);

  // compute C_2
  size_t C_2;
  element_t g_1_r;
  element_init_G1(g_1_r, pairing_);
  element_pow_zn(g_1_r, g_1, r);
  unsigned char *data = (unsigned char *)malloc(element_length_in_bytes(g_1_r));
  element_to_bytes(data, g_1_r);
  C_2 = hash((char *)data) ^ mes;

  /*
   *decrypt
   */
  size_t mes_;
  element_t C_1_alpha;
  element_init_G1(C_1_alpha, pairing_);
  element_pow_zn(C_1_alpha, C_1, alpha);
  element_to_bytes(data, C_1_alpha);
  mes_ = C_2 ^ hash((char *)data);

  if (mes_ == mes) {
    printf("解密结果：%ld\n", mes_);
    printf("解密成功！\n");
  } else {
    printf("解密结果：%ld\n", mes_);
    printf("解密失败！\n");
  }
}

#endif