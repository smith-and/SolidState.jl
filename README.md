# SolidState

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://smith-and.gitlab.io/SolidState/dev)
[![Build Status](https://gitlab.com/smith-and/SolidState/badges/master/pipeline.svg)](https://gitlab.com/smith-and/SolidState/pipelines)
[![Coverage](https://gitlab.com/smith-and/SolidState/badges/master/coverage.svg)](https://gitlab.com/smith-and/SolidState/commits/master)

## Refactoring

I want to change the KinematicDensity to be an abstract type, then have different ChartType specify a subtype of OperatorDensity such that it computes the appropriate operators for the calculation

```julia
abstract type OperatorDensity end
```

Then we would

```julia
struct DataMap{ChartType <: DataChart, KType <: OperatorDensity, LType <: AbstractArray}
    d::Int
    Λ::LType
    K::KType
    chart::ChartType
end
```
