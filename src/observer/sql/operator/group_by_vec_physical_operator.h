/* Copyright (c) 2021 OceanBase and/or its affiliates. All rights reserved.
miniob is licensed under Mulan PSL v2.
You can use this software according to the terms and conditions of the Mulan PSL v2.
You may obtain a copy of Mulan PSL v2 at:
         http://license.coscl.org.cn/MulanPSL2
THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
See the Mulan PSL v2 for more details. */

#pragma once

#include "sql/expr/aggregate_hash_table.h"
#include "sql/operator/physical_operator.h"

/**
 * @brief Group By 物理算子(vectorized)
 * @ingroup PhysicalOperator
 */
class GroupByVecPhysicalOperator : public PhysicalOperator
{
public:
  GroupByVecPhysicalOperator(
      std::vector<std::unique_ptr<Expression>> &&group_by_exprs, std::vector<Expression *> &&expressions);

  virtual ~GroupByVecPhysicalOperator() = default;

  PhysicalOperatorType type() const override { return PhysicalOperatorType::GROUP_BY_VEC; }

  RC open(Trx *trx) override
  {
      ASSERT(children_.size() == 1, "group by operator only support one child, but got %d", children_.size());
      PhysicalOperator &child = *children_[0];
      RC                rc    = child.open(trx);
      if (OB_FAIL(rc)) {
          LOG_INFO("failed to open child operator. rc=%s", strrc(rc));
          return rc;
      }

      while (OB_SUCC(rc = child.next(chunk_))) {

          Chunk group_by_chunk;
          for (size_t group_by_idx = 0; group_by_idx < group_by_expressions_.size(); group_by_idx++) {
              auto column = std::make_unique<Column>();
              group_by_expressions_[group_by_idx]->get_column(chunk_, *column);
              group_by_chunk.add_column(std::move(column), group_by_idx);
          }

          Chunk aggr_chunk;
          for (size_t aggr_idx = 0; aggr_idx < aggregate_expressions_.size(); aggr_idx++) {
              auto column = std::make_unique<Column>();
              value_expressions_[aggr_idx]->get_column(chunk_, *column);
              aggr_chunk.add_column(std::move(column), aggr_idx);
          }

          hash_table_.add_chunk(group_by_chunk, aggr_chunk);
      }
      
      scanner_ = new StandardAggregateHashTable::Scanner(&hash_table_);
      scanner_->open_scan();

      if (rc == RC::RECORD_EOF) {
          rc = RC::SUCCESS;
      }

      return rc;
  }

  RC next(Chunk &chunk) override
  {
      RC rc = scanner_->next(output_chunk_);

      if (OB_FAIL(rc)) {
        return rc;
      }
      rc = chunk.reference(output_chunk_);
      return rc;
  }
  
  RC close() override
  {
    children_[0]->close();
    scanner_->close_scan();
    LOG_INFO("close group by operator");
    return RC::SUCCESS;
  }

private:
  std::vector<std::unique_ptr<Expression>>      group_by_expressions_;      // 分组表达式
  std::vector<Expression *>                     aggregate_expressions_;     // 聚合表达式
  std::vector<Expression *>                     value_expressions_;
  StandardAggregateHashTable                    hash_table_;                // 哈希表
  StandardAggregateHashTable::Scanner           *scanner_;
  Chunk                                         chunk_;                                                       
  Chunk                                         output_chunk_;           
  // bool consumed_;
};