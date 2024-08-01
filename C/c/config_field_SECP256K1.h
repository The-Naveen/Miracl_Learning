/*
 * Copyright (c) 2012-2020 MIRACL UK Ltd.
 *
 * This file is part of MIRACL Core
 * (see https://github.com/miracl/core).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file config_field.h
 * @author Mike Scott
 * @brief Config Curve  Header File
 *
 */

#ifndef CONFIG_FIELD_SECP256K1_H
#define CONFIG_FIELD_SECP256K1_H

#include"core.h"
#include "config_big_256_56.h"

// FP stuff

#define MBITS_SECP256K1 256	/**< Modulus bits */
#define PM1D2_SECP256K1 1     /**< Largest m such that 2^m|(p-1) */
#define MODTYPE_SECP256K1 NOT_SPECIAL  /**< Modulus type */
#define MAXXES_SECP256K1 24     /**< Maximum excess for lazy reduction */
#define QNRI_SECP256K1 0       /**< Small Quadratic Non-Residue */
#define RIADZ_SECP256K1 -11      /**< Z for hash to Curve */
#define RIADZG2A_SECP256K1 0     /**< real part of Z in G2 for Hash to Curve */
#define RIADZG2B_SECP256K1 0     /**< imaginary part of Z in G2 for Hash to Curve */
#define TOWER_SECP256K1 NEGATOWER          /**< Postive or Negative towering */

#endif
