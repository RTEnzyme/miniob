/* Copyright (c) 2021 OceanBase and/or its affiliates. All rights reserved.
miniob is licensed under Mulan PSL v2.
You can use this software according to the terms and conditions of the Mulan PSL v2.
You may obtain a copy of Mulan PSL v2 at:
         http://license.coscl.org.cn/MulanPSL2
THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
See the Mulan PSL v2 for more details. */

#include "sql/operator/group_by_vec_physical_operator.h"
#include <algorithm>

using namespace std;
using namespace common;

GroupByVecPhysicalOperator::GroupByVecPhysicalOperator(
    std::vector<std::unique_ptr<Expression>> &&group_by_exprs, std::vector<Expression *> &&expressions):
    group_by_expressions_(std::move(group_by_exprs)),
    aggregate_expressions_(std::move(expressions)),
    hash_table_(aggregate_expressions_){

    value_expressions_.reserve(aggregate_expressions_.size());

    ranges::for_each(aggregate_expressions_, [this](Expression *expr) {
        auto *      aggregate_expr = static_cast<AggregateExpr *>(expr);
        Expression *child_expr     = aggregate_expr->child().get();
        ASSERT(child_expr != nullptr, "aggregation expression must have a child expression");
        value_expressions_.emplace_back(child_expr);
    });

    // 添加 group by 列
    for (int i = 0; i< static_cast<int>(group_by_expressions_.size()); i++){
        auto group_by_expr = group_by_expressions_[i].get();
        if (group_by_expr->value_type() == AttrType::INTS) {
        output_chunk_.add_column(make_unique<Column>(AttrType::INTS, 4), i);
        } else if (group_by_expr->value_type() == AttrType::FLOATS) {
        output_chunk_.add_column(make_unique<Column>(AttrType::FLOATS, 4), i);
        } else if (group_by_expr->value_type() == AttrType::CHARS) {
        output_chunk_.add_column(make_unique<Column>(AttrType::CHARS, group_by_expr->value_length()), i);
        } else if (group_by_expr->value_type() == AttrType::BOOLEANS) {
        output_chunk_.add_column(make_unique<Column>(AttrType::BOOLEANS, 1), i);
        }
        else {
        ASSERT(false, "not supported groupby type");
        }
    }
    int group_by_expr_size = static_cast<int>(group_by_expressions_.size());

    // 添加 aggr 列
    for (int i = 0; i< static_cast<int>(aggregate_expressions_.size()); i++){
        auto aggregate_expr = static_cast<AggregateExpr *>(aggregate_expressions_[i]);
        if (aggregate_expr->value_type() == AttrType::INTS) {
        output_chunk_.add_column(make_unique<Column>(AttrType::INTS, 4), i+group_by_expr_size);
        } else if (aggregate_expr->value_type() == AttrType::FLOATS) {
        output_chunk_.add_column(make_unique<Column>(AttrType::FLOATS, 4), i+group_by_expr_size);
        }else {
        ASSERT(false, "not supported aggregation type");
        }
    }

};