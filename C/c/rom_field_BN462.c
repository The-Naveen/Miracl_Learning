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
#include "arch.h"
#include "fp_BN462.h"

/* Curve BN462 - Pairing friendly BN curve */

#if CHUNK==16

#error Not supported

#endif

#if CHUNK==32
// Base Bits= 28
const BIG_464_28 Modulus_BN462= {0x138013,0x1B0084,0x24,0xF640000,0xC6FF687,0xF6FF66F,0xFFFFFFF,0xC8020FF,0x2908F41,0xD81,0xA000000,0xB7D9BFC,0x6FF0CF6,0xFFFFFFF,0x23FFF,0x8036012,0x2404};
const BIG_464_28 ROI_BN462= {0x138012,0x1B0084,0x24,0xF640000,0xC6FF687,0xF6FF66F,0xFFFFFFF,0xC8020FF,0x2908F41,0xD81,0xA000000,0xB7D9BFC,0x6FF0CF6,0xFFFFFFF,0x23FFF,0x8036012,0x2404};
const BIG_464_28 R2modp_BN462= {0x88F9612,0xC8B9999,0x247088C,0xDB3ACD5,0xCA792EF,0x6E92E73,0x34F5999,0x4273E13,0x6714A6A,0xCBFE239,0xA0E2617,0x2487CB3,0xBE0EA3C,0x97791E4,0x8A07DE5,0x56CFA97,0x373};
const chunk MConst_BN462= 0x11BB5E5;
const BIG_464_28 CRu_BN462= {0x1A401A,0x1F80A8,0x24,0xF2E0000,0xBBBF3DB,0xF6FF597,0xFFFFFFF,0x22029FF,0x3689F02,0xD81,0x1000000,0x4BD65FC,0x6FF0876,0xFFFFFFF,0x23FFF,0x8036012,0x2404};
const BIG_464_28 Fra_BN462= {0x2575D1A,0xE4BE3FF,0x659DBDE,0xFC7D89,0x93FA118,0xD45D1D,0xCC78D9,0x6217331,0xD547C05,0xC792504,0x9A87E11,0x92ED03A,0x1727085,0xB5A8CC1,0xB40BCFD,0xF4348CB,0x16F};
const BIG_464_28 Frb_BN462= {0x55FF85B,0x204AE09,0x57BEB62,0x9023886,0xD630A20,0xF94F4B,0x84FF0D0,0x14A1A7A,0xB1DBADB,0xB00D8E4,0x1E85F7,0x7493CA0,0x68325B2,0x2008E6B,0xE90EA88,0x7F80940,0x23C3};
const BIG_464_28 SQRTm3_BN462= {0x210022,0x2400CC,0x24,0xEF80000,0xB07F12F,0xF6FF4BF,0xFFFFFFF,0x7C032FF,0x440AEC2,0xD81,0x8000000,0xDFD2FFB,0x6FF03F5,0xFFFFFFF,0x23FFF,0x8036012,0x2404};
#endif

#if CHUNK==64
// Base Bits= 60
const BIG_464_60 Modulus_BN462= {0x401B00840138013L,0x87F640000000002L,0xFFFF6FF66FC6FF6L,0x8F41C8020FFFFFFL,0xD81290L,0xFF0CF6B7D9BFCA0L,0x23FFFFFFFFFF6L,0x24048036012L};
const BIG_464_60 ROI_BN462= {0x401B00840138012L,0x87F640000000002L,0xFFFF6FF66FC6FF6L,0x8F41C8020FFFFFFL,0xD81290L,0xFF0CF6B7D9BFCA0L,0x23FFFFFFFFFF6L,0x24048036012L};
const BIG_464_60 R2modp_BN462= {0x89118D28DC21038L,0x1C24CD524708896L,0x96F6AF594FD13D3L,0xFC17B3AFB34F599L,0x617CBFE0F54B3BCL,0x105034B613F1E2L,0x47E597791E4CB9L,0x12EACA995DAL};
const chunk MConst_BN462= 0x718CE9E711BB5E5L;
const BIG_464_60 CRu_BN462= {0x401F80A801A401AL,0xDBF2E0000000002L,0xFFFF6FF597BBBF3L,0x9F0222029FFFFFFL,0xD81368L,0xFF08764BD65FC10L,0x23FFFFFFFFFF6L,0x24048036012L};
const BIG_464_60 Fra_BN462= {0xEE4BE3FF2575D1AL,0x180FC7D89659DBDL,0x8D90D45D1D93FA1L,0x7C0562173310CC7L,0x87E11C792504D54L,0x72708592ED03A9AL,0xB40BCFDB5A8CC11L,0x16FF4348CBL};
const BIG_464_60 Frb_BN462= {0x2204AE0955FF85BL,0x20902388657BEB6L,0xD00F94F4BD630AL,0xBADB14A1A7A84FFL,0xE85F7B00D8E4B1DL,0x8325B27493CA001L,0xE90EA882008E6B6L,0x23C37F80940L};
const BIG_464_60 SQRTm3_BN462= {0x402400CC0210022L,0x2FEF80000000002L,0xFFFF6FF4BFB07F1L,0xAEC27C032FFFFFFL,0xD81440L,0xFF03F5DFD2FFB80L,0x23FFFFFFFFFF6L,0x24048036012L};
#endif