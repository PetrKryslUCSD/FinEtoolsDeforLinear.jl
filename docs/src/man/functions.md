# Functions

## FEM machines

### Linear deformation

#### Model reduction types

```@autodocs
Modules = [FinEtools, FinEtools.DeforModelRedModule]
Private = true
Order = [:function]
```

#### Base functionality

```@autodocs
Modules = [FinEtools, FinEtools.FEMMDeforLinearBaseModule, FinEtools.FEMMDeforLinearModule, FinEtools.FEMMDeforWinklerModule, FinEtools.FEMMDeforLinearMSModule, FinEtools.FEMMDeforSurfaceDampingModule, FinEtools.FEMMDeforLinearNICEModule, FinEtools.FEMMDeforLinearESNICEModule]
Private = true
Order = [:function]
```

#### Simple FE models

```@autodocs
Modules = [FinEtools, FinEtools.FEMMDeforLinearModule, FinEtools.FEMMDeforWinklerModule,  FinEtools.FEMMDeforSurfaceDampingModule, ]
Private = true
Order = [:function]
```

#### Advanced FE models

```@autodocs
Modules = [FinEtools, FinEtools.FEMMDeforLinearMSModule, FinEtools.FEMMDeforLinearNICEModule, FinEtools.FEMMDeforLinearESNICEModule]
Private = true
Order = [:function]
```

## Algorithms

### Linear deformation

```@autodocs
Modules = [FinEtools, FinEtools.AlgoDeforLinearModule]
Private = true
Order = [:function]
```

## Material models

### Material for deformation, base functionality

```@autodocs
Modules = [FinEtools, FinEtools.MatDeforModule]
Private = true
Order = [:function]
```

### Material models for elasticity

```@autodocs
Modules = [FinEtools, FinEtools.MatDeforLinearElasticModule, FinEtools.MatDeforElastIsoModule, FinEtools.MatDeforElastOrthoModule,]
Private = true
Order = [:function]
```
