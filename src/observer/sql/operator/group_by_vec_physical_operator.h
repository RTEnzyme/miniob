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
      std::vector<std::unique_ptr<Expression>> &&group_by_exprs, std::vector<Expression *> &&expressions){
        aggregate_expressions_ = std::move(expressions);
        group_by_expressions_  = std::move(group_by_exprs);
        // value_expressions_.reserve(aggregate_expressions_.size());

        // ranges::for_each(aggregate_expressions_, [this](Expression *expr) {
        //   auto *      aggregate_expr = static_cast<AggregateExpr *>(expr);
        //   Expression *child_expr     = aggregate_expr->child().get();
        //   ASSERT(child_expr != nullptr, "aggregation expression must have a child expression");
        //   value_expressions_.emplace_back(child_expr);
        // });

        hash_table_ = std::make_unique<StandardAggregateHashTable>(aggregate_expressions_);

        // 添加 group by 列
        for (int i = 0; i< static_cast<int>(group_by_expressions_.size()); i++){
          auto group_by_expr = group_by_expressions_[i].get();
          if (group_by_expr->value_type() == AttrType::INTS) {
            output_chunk_.add_column(make_unique<Column>(group_by_expr->type(), sizeof(int)), i);
          } else if (group_by_expr->value_type() == AttrType::FLOATS) {
            output_chunk_.add_column(make_unique<Column>(group_by_expr->type(), sizeof(float)), i);
          }else {
            ASSERT(false, "not supported aggregation type");
          }
        }

        // 添加聚合列
        for (int i = 0; i< static_cast<int>(aggregate_expressions_.size()); i++){
          auto aggregate_expr = static_cast<AggregateExpr *>(aggregate_expressions_[i]);
          if (aggregate_expr->value_type() == AttrType::INTS) {
            output_chunk_.add_column(make_unique<Column>(aggregate_expressions_[i]->type(), sizeof(int)), i);
          } else if (aggregate_expr->value_type() == AttrType::FLOATS) {
            output_chunk_.add_column(make_unique<Column>(aggregate_expressions_[i]->type(), sizeof(float)), i);
          }else {
            ASSERT(false, "not supported aggregation type");
          }
        }
      };

  virtual ~GroupByVecPhysicalOperator() = default;

  PhysicalOperatorType type() const override { return PhysicalOperatorType::GROUP_BY_VEC; }

  RC open(Trx *trx) override { return RC::UNIMPLENMENT; }
  RC next(Chunk &chunk) override { return RC::UNIMPLENMENT; }
  RC close() override { return RC::UNIMPLENMENT; }

private:
  std::vector<Expression *> aggregate_expressions_;                  // 聚合表达式
  // std::vector<Expression *> value_expressions_;
  std::vector<std::unique_ptr<Expression>> group_by_expressions_;     // 分组表达式
  std::unique_ptr<StandardAggregateHashTable> hash_table_;             //  哈希表
  Chunk chunk_;                                                       
  Chunk output_chunk_;                                                  // 输出 chunk
};