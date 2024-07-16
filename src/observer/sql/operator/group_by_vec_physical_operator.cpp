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

using namespace std;

RC GroupByVecPhysicalOperator::open(Trx *trx)
{
    ASSERT(children_.size() == 1, "group by operator only support one child, but got %d", children_.size());
    PhysicalOperator &child = *children_[0];
    RC                rc    = child.open(trx);
    if (OB_FAIL(rc)) {
        LOG_INFO("failed to open child operator. rc=%s", strrc(rc));
        return rc;
    }

    while (OB_SUCC(rc = child.next(chunk_))) {
        
        Chunk aggr_chunk;
        for (size_t aggr_idx = 0; aggr_idx < aggregate_expressions_.size(); aggr_idx++) {
            auto column = std::make_unique<Column>();
            aggregate_expressions_[aggr_idx]->get_column(chunk_, *column);
            aggr_chunk.add_column(std::move(column), aggr_idx);
        }
                
        Chunk group_by_chunk;
        for (size_t group_by_idx = 0; group_by_idx < group_by_expressions_.size(); group_by_idx++) {
            auto column = std::make_unique<Column>();
            group_by_expressions_[group_by_idx]->get_column(chunk_, *column);
            group_by_chunk.add_column(std::move(column), group_by_idx);
        }
        
        hash_table_->add_chunk(group_by_chunk, aggr_chunk);
    }

    if (rc == RC::RECORD_EOF) {
        rc = RC::SUCCESS;
    }

    return rc;
}

RC GroupByVecPhysicalOperator::next(Chunk &chunk)
{
    if (consumed_) {
        return RC::INTERNAL;
    }
    StandardAggregateHashTable::Scanner scanner(hash_table_.get());
    scanner.open_scan();
    auto rc = scanner.next(output_chunk_);
    rc = chunk.reference(output_chunk_);
    consumed_ = true;
    return rc;
}

RC GroupByVecPhysicalOperator::close()
{
  children_[0]->close();
  LOG_INFO("close group by operator");
  return RC::SUCCESS;
}