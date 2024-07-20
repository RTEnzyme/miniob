/* Copyright (c) 2021 OceanBase and/or its affiliates. All rights reserved.
miniob is licensed under Mulan PSL v2.
You can use this software according to the terms and conditions of the Mulan PSL v2.
You may obtain a copy of Mulan PSL v2 at:
         http://license.coscl.org.cn/MulanPSL2
THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
See the Mulan PSL v2 for more details. */

#include "sql/expr/aggregate_hash_table.h"

// ----------------------------------StandardAggregateHashTable------------------

RC StandardAggregateHashTable::add_chunk(Chunk &groups_chunk, Chunk &aggrs_chunk)
{
    for (int row_idx = 0; row_idx < groups_chunk.rows(); row_idx++) {
      std::vector<Value> group_by_row_value;
      for (int i = 0; i < groups_chunk.column_num(); i++) {
        group_by_row_value.push_back(groups_chunk.get_value(i, row_idx));
      }
      std::vector<Value> aggr_row_value;
      for (int i = 0; i < aggrs_chunk.column_num(); i++) {
        aggr_row_value.push_back(aggrs_chunk.get_value(i, row_idx));
      }
      // aggr_values_ : 哈希表
      auto it = aggr_values_.find(group_by_row_value);
      if (it == end()) {
        // 新建一个 group_by_row_value
        aggr_values_[group_by_row_value] = aggr_row_value;
      } else {
        auto &aggrs = it->second;
        for (int i = 0; i < static_cast<int>(aggrs.size()); i++) {
          auto &aggr = aggrs[i];
          switch (aggrs[i].attr_type())
          {
          case AttrType::INTS:
            aggr.set_int(aggr.get_int() + aggr_row_value[i].get_int());
            break;
          case AttrType::FLOATS:
            aggr.set_float(aggr.get_float() + aggr_row_value[i].get_float());
            break;
          default:
            break;
          }
        }
      }
    }
    return RC::SUCCESS;
}

void StandardAggregateHashTable::Scanner::open_scan()
{
  it_  = static_cast<StandardAggregateHashTable *>(hash_table_)->begin();
  end_ = static_cast<StandardAggregateHashTable *>(hash_table_)->end();
}

RC StandardAggregateHashTable::Scanner::next(Chunk &output_chunk)
{
  if (it_ == end_) {
    return RC::RECORD_EOF;
  }
  while (it_ != end_ && output_chunk.rows() <= output_chunk.capacity()) {
    auto &group_by_values = it_->first;   // std::vector<Value>
    auto &aggrs           = it_->second;  // std::vector<Value>
    for (int i = 0; i < output_chunk.column_num(); i++) {
      auto col_idx = output_chunk.column_ids(i);
      if (col_idx >= static_cast<int>(group_by_values.size())) {
        output_chunk.column(i).append_one((char *)aggrs[col_idx - group_by_values.size()].data());
      } else {
        output_chunk.column(i).append_one((char *)group_by_values[col_idx].data());
      }
    }
    it_++;
  }
  if (it_ == end_) {
    return RC::SUCCESS;
  }

  return RC::SUCCESS;
}

size_t StandardAggregateHashTable::VectorHash::operator()(const vector<Value> &vec) const
{
  size_t hash = 0;
  for (const auto &elem : vec) {
    hash ^= std::hash<string>()(elem.to_string());
  }
  return hash;
}

bool StandardAggregateHashTable::VectorEqual::operator()(const vector<Value> &lhs, const vector<Value> &rhs) const
{
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (size_t i = 0; i < lhs.size(); ++i) {
    if (rhs[i].compare(lhs[i]) != 0) {
      return false;
    }
  }
  return true;
}

// ----------------------------------LinearProbingAggregateHashTable------------------
#ifdef USE_SIMD
template <typename V>
RC LinearProbingAggregateHashTable<V>::add_chunk(Chunk &group_chunk, Chunk &aggr_chunk)
{
  if (group_chunk.column_num() != 1 || aggr_chunk.column_num() != 1) {
    LOG_WARN("group_chunk and aggr_chunk size must be 1.");
    return RC::INVALID_ARGUMENT;
  }
  if (group_chunk.rows() != aggr_chunk.rows()) {
    LOG_WARN("group_chunk and aggr _chunk rows must be equal.");
    return RC::INVALID_ARGUMENT;
  }
  add_batch((int *)group_chunk.column(0).data(), (V *)aggr_chunk.column(0).data(), group_chunk.rows());
  return RC::SUCCESS;
}

template <typename V>
void LinearProbingAggregateHashTable<V>::Scanner::open_scan()
{
  capacity_   = static_cast<LinearProbingAggregateHashTable *>(hash_table_)->capacity();
  size_       = static_cast<LinearProbingAggregateHashTable *>(hash_table_)->size();
  scan_pos_   = 0;
  scan_count_ = 0;
}

template <typename V>
RC LinearProbingAggregateHashTable<V>::Scanner::next(Chunk &output_chunk)
{
  if (scan_pos_ >= capacity_ || scan_count_ >= size_) {
    std::cout << "scan_pos_: " << scan_pos_ << ", capacity_: " << capacity_ << ", scan_count_: " << scan_count_ << ", size_: " << size_ << std::endl;
    return RC::RECORD_EOF;
  }
  auto linear_probing_hash_table = static_cast<LinearProbingAggregateHashTable *>(hash_table_);
  while (scan_pos_ < capacity_ && scan_count_ < size_ && output_chunk.rows() <= output_chunk.capacity()) {
    int key;
    V   value;
    RC  rc = linear_probing_hash_table->iter_get(scan_pos_, key, value);
    if (rc == RC::SUCCESS) {
      output_chunk.column(0).append_one((char *)&key);
      output_chunk.column(1).append_one((char *)&value);
      scan_count_++;
    }
    scan_pos_++;
  }
  return RC::SUCCESS;
}

template <typename V>
void LinearProbingAggregateHashTable<V>::Scanner::close_scan()
{
  capacity_   = -1;
  size_       = -1;
  scan_pos_   = -1;
  scan_count_ = 0;
}

template <typename V>
RC LinearProbingAggregateHashTable<V>::get(int key, V &value)
{
  RC  rc          = RC::SUCCESS;
  int index       = (key % capacity_ + capacity_) % capacity_;
  int iterate_cnt = 0;
  while (true) {
    if (keys_[index] == EMPTY_KEY) {
      rc = RC::NOT_EXIST;
      break;
    } else if (keys_[index] == key) {
      value = values_[index];
      break;
    } else {
      index += 1;
      index %= capacity_;
      iterate_cnt++;
      if (iterate_cnt > capacity_) {
        rc = RC::NOT_EXIST;
        break;
      }
    }
  }
  return rc;
}

template <typename V>
RC LinearProbingAggregateHashTable<V>::iter_get(int pos, int &key, V &value)
{
  RC rc = RC::SUCCESS;
  if (keys_[pos] == LinearProbingAggregateHashTable<V>::EMPTY_KEY) {
    rc = RC::NOT_EXIST;
  } else {
    key   = keys_[pos];
    value = values_[pos];
  }
  return rc;
}

template <typename V>
void LinearProbingAggregateHashTable<V>::aggregate(V *value, V value_to_aggregate)
{
  if (aggregate_type_ == AggregateExpr::Type::SUM) {
    *value += value_to_aggregate;
  } else {
    ASSERT(false, "unsupported aggregate type");
  }
}

template <typename V>
void LinearProbingAggregateHashTable<V>::resize()
{
  capacity_ *= 2;
  std::vector<int> new_keys(capacity_);
  std::vector<V>   new_values(capacity_);

  for (size_t i = 0; i < keys_.size(); i++) {
    auto &key   = keys_[i];
    auto &value = values_[i];
    if (key != EMPTY_KEY) {
      int index = (key % capacity_ + capacity_) % capacity_;
      while (new_keys[index] != EMPTY_KEY) {
        index = (index + 1) % capacity_;
      }
      new_keys[index]   = key;
      new_values[index] = value;
    }
  }

  keys_   = std::move(new_keys);
  values_ = std::move(new_values);
}

template <typename V>
void LinearProbingAggregateHashTable<V>::resize_if_need()
{
  if (size_ >= capacity_ / 2) {
    resize();
  }
}

__m256i calculate_hash_and_mod(__m256i vec_keys, int table_size) {
  // (key % capacity_ + capacity_) % capacity_;
    // 计算 table_size - 1
    int mask = table_size - 1;
    __m256i vec_mask = _mm256_set1_epi32(mask);
    __m256i vec_table_size = _mm256_set1_epi32(table_size);

    // key & (table_size - 1)
    __m256i vec_mod = _mm256_and_si256(vec_keys, vec_mask);
    
    // 防止负值出现，如果出现则加上 table_size
    __m256i vec_adjusted = _mm256_add_epi32(vec_mod, vec_table_size);
    __m256i vec_result = _mm256_and_si256(vec_adjusted, vec_mask);
    
    return vec_result;
}

template <typename V>
void scatter_epi32(V *keys_ptr, __m256i vec_idx, __m256i vec_keys_get) {
    int *idx = (int *)&vec_idx;
    int *keys = (int *)&vec_keys_get;

    // 按索引将值存储到 keys_ptr 中
    for (int i = 0; i < SIMD_WIDTH; i++) {
        keys_ptr[idx[i]] += keys[i];
    }
}

__m256i compare_vectors(__m256i vec_1, __m256i vec_2) {
    // 使用比较指令比较两个向量的元素
    __m256i vec_cmp = _mm256_cmpeq_epi32(vec_1, vec_2);

    return vec_cmp;
}

void print__mm256i(__m256i vec, string prefix) {
    int *ptr = (int *)&vec;
    std::cout << prefix << ": ";
    for (int i = 0; i < SIMD_WIDTH; i++) {
        std::cout << ptr[i] << " ";
    }
    std::cout << std::endl;
}


template <typename V>
void LinearProbingAggregateHashTable<V>::add_batch(int *input_keys, V *input_values, int len)
{
  // TODO
  __m256i inv = _mm256_set1_epi32(-1);
  __m256i off = _mm256_set1_epi32(0);
  int i = 0;
  int tmp_keys[SIMD_WIDTH];
  V tmp_values[SIMD_WIDTH];

  for (; i + SIMD_WIDTH <= len;) {
    // std::cout << "i: " << i << std::endl;
    // 读数据
    selective_load(input_keys, i, tmp_keys, inv);
    __m256i vec_input_keys = _mm256_loadu_si256((__m256i const *)tmp_keys);
    
    selective_load(input_values, i, tmp_values, inv);
    __m256i vec_input_values = _mm256_loadu_si256((__m256i const *)tmp_values);



    // 计算哈希
    __m256i vec_idx = calculate_hash_and_mod(vec_input_keys, capacity_);
    // vec_idx = vec_idx + off
    vec_idx = _mm256_add_epi32(vec_idx, off);

    // 从 keys_ 和 values_ 中取出 key 和 value
    int* keys_ptr = keys_.data();
    __m256i vec_keys_get = _mm256_i32gather_epi32(keys_ptr, vec_idx, 4);


    // V* values_ptr = values_.data();
    // __m256i vec_values_get = _mm256_i32gather_epi32(values_ptr, vec_idx, 4); // float 和 int 都是 4 字节

    // 检测是否为 EMPTY_KEY => -1
    __m256i vec_empty_keys = _mm256_set1_epi32(EMPTY_KEY);
    __m256i vec_cmp_1 = compare_vectors(vec_keys_get, vec_empty_keys);


    // 统计 vec_cmp_1 中的非零元素个数
    const int *cmp_1_ptr = reinterpret_cast<const int *>(&vec_cmp_1);
    size_ -= mm256_sum_epi32(cmp_1_ptr, SIMD_WIDTH);

    // 检测两个 key 是否相同 => -1
    __m256i vec_cmp_2 = compare_vectors(vec_keys_get, vec_input_keys);


    // vec_cmp_1 | vec_cmp_2 => -1 / 0
    __m256i vec_cmp = _mm256_or_si256(vec_cmp_1, vec_cmp_2);


    // TODO
    // int *keys_get_ptr = reinterpret_cast<int *>(&vec_keys_get); // 哈希表中的 key
    int *idx_ptr = reinterpret_cast<int *>(&vec_idx); // 哈希表中的 idx
    __m256i vec_cmp_3 = _mm256_setzero_si256();
    int *cmp_ptr_3 = reinterpret_cast<int *>(&vec_cmp_3);
    for (int j = 0; j < SIMD_WIDTH; j++) {
      if (keys_[idx_ptr[j]] == EMPTY_KEY) {
        keys_[idx_ptr[j]] = tmp_keys[j];
        cmp_ptr_3[j] = -1;
      } else if (keys_[idx_ptr[j]] == tmp_keys[j]){
        cmp_ptr_3[j] = -1;
      } else {
        cmp_ptr_3[j] = 0;
      }
      // std::cout << "idx: " << idx_ptr[j] << ", key: " << keys_[idx_ptr[j]] << ", value: " << values_[idx_ptr[j]] << endl;
    }

    // 计算 vec_cmp 
    vec_cmp = _mm256_and_si256(vec_cmp, vec_cmp_3);


    inv = vec_cmp;


    // 将 inv 转换为 const int *
    const int *inv_ptr = reinterpret_cast<const int *>(&inv);
    i -= mm256_sum_epi32(inv_ptr, SIMD_WIDTH);

    // std::cout << "i2: " << i << std::endl;

    // off + (1+vec_cmp)
    // 清空 inv 为 -1 的位
    off = _mm256_blendv_epi8(off, _mm256_set1_epi32(0), vec_cmp);
    off = _mm256_add_epi32(off, _mm256_add_epi32(_mm256_set1_epi32(1), vec_cmp));

    // 对 vec_cmp 为 1 的元素，则 vec_aggregate 相应元素 +vec_input_values；否则vec_aggregate 相应元素 +0
    __m256i vec_zero = _mm256_setzero_si256();
    __m256i vec_aggregate = _mm256_blendv_epi8(vec_zero, vec_input_values, vec_cmp);

    scatter_epi32(values_.data(), vec_idx, vec_aggregate);
  }
  // int *inv_ptr = reinterpret_cast<int *>(&inv);
  // for (int j = 0; j < SIMD_WIDTH; j++) {
  //   if (inv_ptr[j] == 0) {
  //     int key = tmp_keys[j];
  //     V value = tmp_values[j];
  //     int index = (key % capacity_ + capacity_) % capacity_;
  //     while (true) {
  //       if (keys_[index] == EMPTY_KEY) {
  //         keys_[index] = key;
  //         values_[index] = value;
  //         size_ +=1;
  //         break;
  //       } else if (keys_[index] == key) {
  //         aggregate(&values_[index], value);
  //         // std::cout << "idx: " << index << ", key: " << keys_[index] << ", value: " << values_[index] << std::endl;
  //         size_ +=1;
  //         break;
  //       } else {
  //         index = (index + 1) % capacity_;
  //       }
  //     }
  //     i += 1;
  //   }   
  // }

  // int left_keys[SIMD_WIDTH];
  // V left_values[SIMD_WIDTH];
  // inv = _mm256_set1_epi32(-1);
  // // std::cout << "i: " << i << std::endl;
  // selective_load_left(input_keys, i, left_keys, inv, len);
  // selective_load_left(input_values, i, left_values, inv, len);
  // for (int j = 0; j < SIMD_WIDTH; j++) {
  //   int key = left_keys[j];
  //   if (key == EMPTY_KEY) {
  //     continue;
  //   }
  //   V value = left_values[j];
  //   int index = (key % capacity_ + capacity_) % capacity_;
  //   while (true) {
  //     if (keys_[index] == EMPTY_KEY) {
  //       keys_[index] = key;
  //       values_[index] = value;
  //       size_ +=1;
  //       break;
  //     } else if (keys_[index] == key) {
  //       aggregate(&values_[index], value);
  //       // std::cout << "idx: " << index << ", key: " << keys_[index] << ", value: " << values_[index] << std::endl;
  //       size_ +=1;
  //       break;
  //     } else {
  //       index = (index + 1) % capacity_;
  //     }
  //   }
  // }
  
  
  // 标量线性探测
  // for (int i=0;i < len; i++) {
  //   int key = input_keys[i];
  //   V value = input_values[i];
  //   int index = (key % capacity_ + capacity_) % capacity_;
  //   while (true) {
  //     if (keys_[index] == EMPTY_KEY) {
  //       keys_[index] = key;
  //       values_[index] = value;
  //       size_ +=1;
  //       break;
  //     } else if (keys_[index] == key) {
  //       aggregate(&values_[index], value);
  //       size_ +=1;
  //       break;
  //     } else {
  //       index = (index + 1) % capacity_;
  //     }
  //   }
  // }

  // size_ = 8;
  // std::cout << "current size_: " << size_ << std::endl;
  

  // inv (invalid) 表示是否有效，inv[i] = -1 表示有效，inv[i] = 0 表示无效。
  // key[SIMD_WIDTH],value[SIMD_WIDTH] 表示当前循环中处理的键值对。
  // off (offset) 表示线性探测冲突时的偏移量，key[i] 每次遇到冲突键，则off[i]++，如果key[i] 已经完成聚合，则off[i] = 0，
  // i = 0 表示selective load 的起始位置。
  // inv 全部初始化为 -1
  // off 全部初始化为 0

  // for (; i + SIMD_WIDTH <= len;) {
    // 1: 根据 `inv` 变量的值，从 `input_keys` 中 `selective load` `SIMD_WIDTH` 个不同的输入键值对。
    // 2. 计算 i += |inv|, `|inv|` 表示 `inv` 中有效的个数 
    // 3. 计算 hash 值，
    // 4. 根据聚合类型（目前只支持 sum），在哈希表中更新聚合结果。如果本次循环，没有找到key[i] 在哈希表中的位置，则不更新聚合结果。
    // 5. gather 操作，根据 hash 值将 keys_ 的 gather 结果写入 table_key 中。
    // 6. 更新 inv 和 off。如果本次循环key[i] 聚合完成，则inv[i]=-1，表示该位置在下次循环中读取新的键值对。
    // 如果本次循环 key[i] 未在哈希表中聚合完成（table_key[i] != key[i]），则inv[i] = 0，表示该位置在下次循环中不需要读取新的键值对。
    // 如果本次循环中，key[i]聚合完成，则off[i] 更新为 0，表示线性探测偏移量为 0，key[i] 未完成聚合，则off[i]++,表示线性探测偏移量加 1。
  // }
  //7. 通过标量线性探测，处理剩余键值对

  // resize_if_need();
}

template <typename V>
const int LinearProbingAggregateHashTable<V>::EMPTY_KEY = 0xffffffff;
template <typename V>
const int LinearProbingAggregateHashTable<V>::DEFAULT_CAPACITY = 16384;

template class LinearProbingAggregateHashTable<int>;
template class LinearProbingAggregateHashTable<float>;
#endif