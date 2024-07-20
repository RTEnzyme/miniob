/* Copyright (c) 2021 OceanBase and/or its affiliates. All rights reserved.
miniob is licensed under Mulan PSL v2.
You can use this software according to the terms and conditions of the Mulan PSL v2.
You may obtain a copy of Mulan PSL v2 at:
         http://license.coscl.org.cn/MulanPSL2
THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
See the Mulan PSL v2 for more details. */

#include <stdint.h>
#include "common/math/simd_util.h"

#if defined(USE_SIMD)

int mm256_extract_epi32_var_indx(const __m256i vec, const unsigned int i)
{
  __m128i idx = _mm_cvtsi32_si128(i);
  __m256i val = _mm256_permutevar8x32_epi32(vec, _mm256_castsi128_si256(idx));
  return _mm_cvtsi128_si32(_mm256_castsi256_si128(val));
}

int mm256_sum_epi32(const int *values, int size)
{
    // TODO
    // 定义用于存储 部分和 的AVX向量
    __m256i vec_sum = _mm256_setzero_si256();

    // 每次处理8个整数
    int i = 0;
    for (; i <= size - 8; i += 8) {
        __m256i vec_values = _mm256_loadu_si256((__m256i const *)(values + i));
        vec_sum = _mm256_add_epi32(vec_sum, vec_values);
    }

    // 将 AVX 向量中的部分和提取到数组中
    alignas(32) int sums[8];
    _mm256_store_si256((__m256i *)sums, vec_sum);

    // 将数组中的部分和累加到最终结果
    int sum = 0;
    for (int j = 0; j < 8; ++j) {
        sum += sums[j];
    }

    // 处理剩余未对齐的整数
    for (; i < size; ++i) {
        sum += values[i];
    }

    return sum;
}

float mm256_sum_ps(const float *values, int size)
{
    // TODO
    // 定义用于存储部分和的AVX向量
    __m256 vec_sum = _mm256_setzero_ps();

    // 每次处理8个浮点数
    int i = 0;
    for (; i <= size - 8; i += 8) {
        __m256 vec_values = _mm256_loadu_ps(values + i);
        vec_sum = _mm256_add_ps(vec_sum, vec_values);
    }

    // 将AVX向量中的部分和提取到数组中
    alignas(32) float sums[8];
    _mm256_store_ps(sums, vec_sum);

    // 将数组中的部分和累加到最终结果
    float sum = 0;
    for (int j = 0; j < 8; ++j) {
        sum += sums[j];
    }

    // 处理剩余未对齐的浮点数
    for (; i < size; ++i) {
        sum += values[i];
    }

    return sum;
}

template <typename V>
void selective_load(V *memory, int offset, V *vec, __m256i &inv)
{
  int *inv_ptr = reinterpret_cast<int *>(&inv);
  for (int i = 0; i < SIMD_WIDTH; i++) {
    if (inv_ptr[i] == -1) {
      vec[i] = memory[offset++];
    } 
  }
}
template <typename V>
void selective_load_left(V *memory, int offset, V *vec, __m256i &inv, int len)
{
  int *inv_ptr = reinterpret_cast<int *>(&inv);
  for (int i = 0; i < SIMD_WIDTH; i++) {
    // 需要读
    if (inv_ptr[i] == -1) {
      if (offset + 1 < len) {
        vec[i] = memory[offset++];
      } else {
        vec[i] = 0xffffffff;
      }
    }
  }
}
template void selective_load<uint32_t>(uint32_t *memory, int offset, uint32_t *vec, __m256i &inv);
template void selective_load<int>(int *memory, int offset, int *vec, __m256i &inv);
template void selective_load<float>(float *memory, int offset, float *vec, __m256i &inv);


template void selective_load_left<uint32_t>(uint32_t *memory, int offset, uint32_t *vec, __m256i &inv, int len);
template void selective_load_left<int>(int *memory, int offset, int *vec, __m256i &inv, int len);
template void selective_load_left<float>(float *memory, int offset, float *vec, __m256i &inv, int len);

#endif