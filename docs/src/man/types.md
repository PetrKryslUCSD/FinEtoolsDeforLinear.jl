# Types

## FEM machines

### Linear deformation

```@autodocs
Modules = [FinEtools, FinEtools.DeforModelRedModule, FinEtools.FEMMDeforLinearBaseModule, FinEtools.FEMMDeforLinearModule, FinEtools.FEMMDeforWinklerModule, FinEtools.FEMMDeforLinearMSModule, FinEtools.FEMMDeforSurfaceDampingModule, FinEtools.FEMMDeforLinearNICEModule, FinEtools.FEMMDeforLinearESNICEModule]
Private = true
Order = [:type]
```

## Material models

### Material for deformation, base functionality

```@autodocs
Modules = [FinEtools, FinEtools.MatDeforModule]
Private = true
Order = [:type]
```

### Material models for elasticity

```@autodocs
Modules = [FinEtools, FinEtools.MatDeforLinearElasticModule, FinEtools.MatDeforElastIsoModule, FinEtools.MatDeforElastOrthoModule,]
Private = true
Order = [:type]
```
