#ifndef CRAMER_SHOUP_H
#define CRAMER_SHOUP_H

#include "curve/curve.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <pbc/pbc.h>
#include <pbc/pbc_field.h>
#include <pbc/pbc_pairing.h>
#include <string>

void test_cramer_shoup() {
  printf("------------------------- Cramer-Shoup -------------------------\n");

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
   * key_gen:
   * sk = (alpha_1, alpha_2, beta_1, beta_2, gamma_1, gamma_2)
   * pk = (u = g_1^alpha_1 * g_2^alpha_2, v = g_1^beta_1 * g_2^beta_2,
   *       h = g_1^gamma_1 * g_2^gamma_2)
   */
  element_t g_1, g_2;
  element_init_G1(g_1, pairing_);
  element_init_G1(g_2, pairing_);
  element_random(g_1);
  element_random(g_2);

  element_t alpha_1, alpha_2, beta_1, beta_2, gamma_1, gamma_2;
  element_init_Zr(alpha_1, pairing_);
  element_init_Zr(alpha_2, pairing_);
  element_init_Zr(beta_1, pairing_);
  element_init_Zr(beta_2, pairing_);
  element_init_Zr(gamma_1, pairing_);
  element_init_Zr(gamma_2, pairing_);
  element_random(alpha_1);
  element_random(alpha_2);
  element_random(beta_1);
  element_random(beta_2);
  element_random(gamma_1);
  element_random(gamma_2);

  element_t u, v, h, tmp;
  element_init_G1(u, pairing_);
  element_init_G1(v, pairing_);
  element_init_G1(h, pairing_);
  element_init_G1(tmp, pairing_);

  // compute u
  element_pow_zn(u, g_1, alpha_1);
  element_pow_zn(tmp, g_2, alpha_2);
  element_mul(u, u, tmp);

  // compute v
  element_pow_zn(v, g_1, beta_1);
  element_pow_zn(tmp, g_2, beta_2);
  element_mul(v, v, tmp);

  // compute h
  element_pow_zn(h, g_1, gamma_1);
  element_pow_zn(tmp, g_2, gamma_2);
  element_mul(h, h, tmp);

  /*
   * encrypt:
   * ct = (C_1 = g1^r, C_2 = g2^r, C_3 = h^r * m, C_4 = u^r * v^{wr})
   */
  element_t mes;
  element_init_G1(mes, pairing_);
  element_random(mes);
  element_printf("明文：%B\n", mes);
  printf("\n");

  element_t r;
  element_init_Zr(r, pairing_);
  element_random(r);

  element_t C_1, C_2, C_3, C_4, w;
  element_init_G1(C_1, pairing_);
  element_init_G1(C_2, pairing_);
  element_init_G1(C_3, pairing_);
  element_init_G1(C_4, pairing_);
  element_init_Zr(w, pairing_);

  // compute C_1
  element_pow_zn(C_1, g_1, r);

  // compute C_2
  element_pow_zn(C_2, g_2, r);

  // compute C_3;
  element_pow_zn(C_3, h, r);
  element_mul(C_3, C_3, mes);

  // compute C_4
  size_t length = element_length_in_bytes(C_1);
  unsigned char *data = (unsigned char *)malloc(length * 3);
  element_to_bytes(data, C_1);
  element_to_bytes(data + length, C_2);
  element_to_bytes(data + 2 * length, C_3);
  element_set_si(w, hash((char *)data));
  element_pow_zn(C_4, u, r);
  element_pow_zn(tmp, v, w);
  element_pow_zn(tmp, tmp, r);
  element_mul(C_4, C_4, tmp);

  /*
   *decrypt
   */
  element_t w_;
  element_init_Zr(w_, pairing_);
  element_to_bytes(data, C_1);
  element_to_bytes(data + length, C_2);
  element_to_bytes(data + 2 * length, C_3);
  element_set_si(w_, hash((char *)data));

  // verify
  element_t C_4_;
  element_init_G1(C_4_, pairing_);
  element_pow_zn(C_4_, C_1, alpha_1);
  element_pow_zn(tmp, C_1, beta_1);
  element_pow_zn(tmp, tmp, w_);
  element_mul(C_4_, C_4_, tmp);
  element_pow_zn(tmp, C_2, alpha_2);
  element_mul(C_4_, C_4_, tmp);
  element_pow_zn(tmp, C_2, beta_2);
  element_pow_zn(tmp, tmp, w_);
  element_mul(C_4_, C_4_, tmp);
  if (!element_cmp(C_4, C_4_)) {
    printf("C_4校验成功！\n");
  } else {
    printf("C_4校验失败！\n");
  }

  // decrypt
  element_t neg_gamma_1, neg_gamma_2;
  element_init_Zr(neg_gamma_1, pairing_);
  element_init_Zr(neg_gamma_2, pairing_);
  element_neg(neg_gamma_1, gamma_1);
  element_neg(neg_gamma_2, gamma_2);

  element_t C_1_neg_gamma_1, C_2_neg_gamma_2;
  element_init_G1(C_1_neg_gamma_1, pairing_);
  element_init_G1(C_2_neg_gamma_2, pairing_);
  element_pow_zn(C_1_neg_gamma_1, C_1, neg_gamma_1);
  element_pow_zn(C_2_neg_gamma_2, C_2, neg_gamma_2);

  element_t mes_;
  element_init_G1(mes_, pairing_);
  element_mul(mes_, C_3, C_1_neg_gamma_1);
  element_mul(mes_, mes_, C_2_neg_gamma_2);

  if (!element_cmp(mes, mes_)) {
    element_printf("解密结果：%B\n", mes_);
    printf("解密成功！\n");
  } else {
    element_printf("解密结果：%B\n", mes_);
    printf("解密失败！\n");
  }
}

#endif