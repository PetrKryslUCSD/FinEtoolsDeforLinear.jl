# Reference Manual

## Simple FE model (volume)

```@autodocs
Modules = [FinEtools, FinEtoolsDeforLinear.FEMMDeforLinearBaseModule, FinEtoolsDeforLinear.FEMMDeforLinearModule,  ]
Private = true
Order = [:function, :type]
```

## Simple FE models (surface)

```@autodocs
Modules = [FinEtools, FinEtoolsDeforLinear.FEMMDeforWinklerModule,  FinEtoolsDeforLinear.FEMMDeforSurfaceDampingModule, ]
Private = true
Order = [:function, :type]
```

## Advanced: Mean-strain FEM

```@autodocs
Modules = [FinEtools, FinEtoolsDeforLinear.FEMMDeforLinearMSModule, ]
Private = true
Order = [:function, :type]
```

## Advanced: Nodal integration

```@autodocs
Modules = [FinEtools, FinEtoolsDeforLinear.FEMMDeforLinearNICEModule, FinEtoolsDeforLinear.FEMMDeforLinearESNICEModule, ]
Private = true
Order = [:function, :type]
```

## Advanced: Incompatible modes

```@autodocs
Modules = [FinEtools, FinEtoolsDeforLinear.FEMMDeforLinearIMModule]
Private = true
Order = [:function, :type]
```

## Algorithms

### Linear deformation

```@autodocs
Modules = [FinEtools, FinEtoolsDeforLinear.AlgoDeforLinearModule]
Private = true
Order = [:function]
```

## Material models

### Material for deformation, base functionality

```@autodocs
Modules = [FinEtools, FinEtoolsDeforLinear.MatDeforModule, ]
Private = true
Order = [:function, :type]
```

### Elasticity

```@autodocs
Modules = [FinEtools, FinEtoolsDeforLinear.MatDeforLinearElasticModule, ]
Private = true
Order = [:function, :type]
```

### Isotropic elasticity

```@autodocs
Modules = [FinEtools, FinEtoolsDeforLinear.MatDeforElastIsoModule, ]
Private = true
Order = [:function, :type]
```

### Orthotropic elasticity

```@autodocs
Modules = [FinEtools,  FinEtoolsDeforLinear.MatDeforElastOrthoModule,]
Private = true
Order = [:function, :type]
```

## Modules

```@docs
FinEtoolsDeforLinear.FEMMDeforLinearESNICEModule
FinEtoolsDeforLinear.FEMMDeforLinearBaseModule
FinEtoolsDeforLinear.FinEtoolsDeforLinear
FinEtoolsDeforLinear.FEMMDeforLinearMSModule
FinEtoolsDeforLinear.MatDeforModule
FinEtoolsDeforLinear.FEMMDeforLinearModule
FinEtoolsDeforLinear.AlgoDeforLinearModule
FinEtoolsDeforLinear.FEMMDeforLinearNICEModule
FinEtoolsDeforLinear.FEMMDeforSurfaceDampingModule
FinEtoolsDeforLinear.MatDeforElastIsoModule
FinEtoolsDeforLinear.MatDeforLinearElasticModule
FinEtoolsDeforLinear.FEMMDeforWinklerModule
FinEtoolsDeforLinear.MatDeforElastOrthoModule
FinEtoolsDeforLinear.FEMMDeforLinearIMModule
```