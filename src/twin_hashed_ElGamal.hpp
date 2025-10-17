#ifndef TWIN_HASHED_ELGAMAL_H
#define TWIN_HASHED_ELGAMAL_H

#include "curve/curve.hpp"
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <pbc/pbc.h>
#include <pbc/pbc_field.h>
#include <pbc/pbc_pairing.h>
#include <string>

void test_twin_hashed_ElGamal() {
  printf("------------------------- twin hashed ElGamal "
         "-------------------------\n");

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
   * key_gen: sample alpha,beta,
   * sk = (alpha,beta), pk = (g_1=g^alpha,g_2=g^beta)
   */
  element_t alpha, beta, g_1, g_2;
  element_init_Zr(alpha, pairing_);
  element_init_Zr(beta, pairing_);
  element_init_G1(g_1, pairing_);
  element_init_G1(g_2, pairing_);
  element_random(alpha);
  element_pow_zn(g_1, g, alpha);
  element_random(beta);
  element_pow_zn(g_2, g, beta);

  /*
   * encrypt:
   * ct = (C_1 = g^r, C_2 = H(g_1^r||g_2^r) xor m)
   */
  size_t mes = 12345;
  printf("明文：%ld\n", mes);

  // sample r
  element_t r;
  element_init_Zr(r, pairing_);
  element_random(r);

  // compute C_1
  element_t C_1;
  element_init_G1(C_1, pairing_);
  element_pow_zn(C_1, g, r);

  // compute C_2
  size_t C_2;
  element_t g_1_r, g_2_r;
  element_init_G1(g_1_r, pairing_);
  element_init_G1(g_2_r, pairing_);
  element_pow_zn(g_1_r, g_1, r);
  element_pow_zn(g_2_r, g_2, r);
  size_t length = element_length_in_bytes(g_1_r);
  unsigned char *data = (unsigned char *)malloc(length * 2);
  element_to_bytes(data, g_1_r);
  element_to_bytes(data + length, g_2_r);
  C_2 = hash((char *)data) ^ mes;

  /*
   *decrypt
   */
  size_t mes_;
  element_t C_1_alpha, C_1_beta;
  element_init_G1(C_1_alpha, pairing_);
  element_init_G1(C_1_beta, pairing_);
  element_pow_zn(C_1_alpha, C_1, alpha);
  element_pow_zn(C_1_beta, C_1, beta);
  element_to_bytes(data, C_1_alpha);
  element_to_bytes(data + length, C_1_beta);
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