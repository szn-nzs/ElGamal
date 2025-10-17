#ifndef ELGAMAL_H
#define ELGAMAL_H

#include "curve/curve.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <pbc/pbc.h>
#include <pbc/pbc_field.h>
#include <pbc/pbc_pairing.h>
#include <string>

void test_ElGamal() {
  printf("------------------------- ElGamal -------------------------\n");

  /*
   * sys_gen: choose a cyclic group and a generator g
   */
  pairing_t pairing_;
  pbc_param_t par;
  pbc_param_init_set_str(par, a_param.c_str());
  pairing_init_pbc_param(pairing_, par);

  element_t g;
  element_init_G1(g, pairing_);
  element_random(g);

  /*
   * key_gen: sample alpha, sk = alpha, pk = g_1=g^alpha
   */
  element_t alpha, g_1;
  element_init_Zr(alpha, pairing_);
  element_init_G1(g_1, pairing_);
  element_random(alpha);
  element_pow_zn(g_1, g, alpha);

  /*
   * encrypt: ct = (C_1 = g^r, C_2 = g_1^r * m)
   */
  element_t mes;
  element_init_G1(mes, pairing_);
  element_random(mes);
  element_printf("明文：%B\n", mes);

  element_t r;
  element_init_Zr(r, pairing_);
  element_random(r);

  // compute C_1
  element_t C_1;
  element_init_G1(C_1, pairing_);
  element_pow_zn(C_1, g, r);

  // compute C_2
  element_t g_1_r;
  element_init_G1(g_1_r, pairing_);
  element_pow_zn(g_1_r, g_1, r);
  element_t C_2;
  element_init_G1(C_2, pairing_);
  element_mul(C_2, g_1_r, mes);

  /*
   *decrypt
   */
  element_t neg_alpha, C_1_neg_alpha;
  element_init_Zr(neg_alpha, pairing_);
  element_init_G1(C_1_neg_alpha, pairing_);
  element_neg(neg_alpha, alpha);
  element_pow_zn(C_1_neg_alpha, C_1, neg_alpha);

  element_t mes_;
  element_init_G1(mes_, pairing_);
  element_mul(mes_, C_2, C_1_neg_alpha);

  if (!element_cmp(mes, mes_)) {
    element_printf("解密结果：%B\n", mes_);
    printf("解密成功！\n");
  } else {
    element_printf("解密结果：%B\n", mes_);
    printf("解密失败！\n");
  }
}

#endif