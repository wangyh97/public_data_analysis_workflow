# Data contracts

A data contract is the interface between processing, analysis and plotting.
It defines the grain (one row represents what), required columns, units,
missing-value rules and unique key. It does **not** demand one universal table.

Contract families:

| Contract | Grain | Typical consumer |
| --- | --- | --- |
| Expression matrix | feature × specimen | bulk correlation/survival |
| Patient observation | patient × timepoint × compartment × feature | single-cell correlations/trajectories |
| Statistical result | one hypothesis | heatmaps, forest plots, annotations |
| Communication result | sample/patient × sender × receiver × pathway | CellChat figures |

See [expression](expression.md), [observations](observations.md), [statistical results](statistical_results.md)
and [communication](communication.md). New exports should record a
`contract_version` and expression unit in their manifest. The CRC example
currently exports v1 patient observations and statistical results; its frozen
legacy outputs predate these contracts. Breaking column or meaning changes
increment the major contract version; adding optional fields increments the
minor version. The v1 CSV validator checks required columns and values; it
does not yet implement version negotiation.
