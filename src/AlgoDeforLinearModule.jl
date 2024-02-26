"""
Module for algorithms used in linear deformation models.
"""
module AlgoDeforLinearModule

__precompile__(true)

using FinEtools.FTypesModule: FDataDict
using FinEtools.AlgoBaseModule: dcheck!
using Arpack: eigs
using SparseArrays: spzeros
using FinEtools.FieldModule:
    AbstractField,
    ndofs,
    setebc!,
    numberdofs!,
    applyebc!,
    scattersysvec!,
    nalldofs,
    nfreedofs,
    gathersysvec
using FinEtools.NodalFieldModule: NodalField, nnodes
import FinEtools.FEMMBaseModule:
    associategeometry!, distribloads, fieldfromintegpoints, elemfieldfromintegpoints
using FinEtoolsDeforLinear.FEMMDeforLinearBaseModule:
    stiffness, mass, thermalstrainloads, inspectintegpoints
using FinEtoolsDeforLinear.FEMMDeforLinearMSModule:
    stiffness, mass, thermalstrainloads, inspectintegpoints
using FinEtools.DeforModelRedModule: stresscomponentmap
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtools.ForceIntensityModule: ForceIntensity
using FinEtools.MeshModificationModule: meshboundary
using FinEtools.MeshExportModule.VTK: vtkexportmesh
using LinearAlgebra:  mul!, norm, eigen, qr, dot, cholesky, Symmetric
my_A_mul_B!(C, A, B) = mul!(C, A, B)

"""
    AlgoDeforLinearModule.linearstatics(modeldata::FDataDict)

Algorithm for static linear deformation (stress) analysis.

# Argument

`modeldata` = dictionary with values for keys

  - `"fens"`  = finite element node set
  - `"regions"`  = array of region dictionaries
  - `"essential_bcs"` = array of essential boundary condition dictionaries
  - `"traction_bcs"` = array of traction boundary condition dictionaries
  - `"temperature_change"` = dictionary of data for temperature change

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:

  - `"femm"` = finite element model machine (mandatory);

For essential boundary conditions (optional) each dictionary
would hold

  - `"displacement"` = fixed (prescribed) displacement (scalar),  or
    a function with signature function `w = f(x)`.
    If this value is not given, zero displacement is assumed.
  - `"component"` = which component is prescribed  (1, 2, 3)?
  - `"node_list"` = list of nodes on the boundary to which the condition applies
    (mandatory)

For traction boundary conditions (optional) each dictionary
would hold key-value pairs

  - `"femm"` = finite element model machine (mandatory);
  - `"traction_vector"` = traction vector,  either  a constant numerical
    vector, or  a function to be used to construct a `ForceIntensity`
    object, or it could be the `ForceIntensity` object itself.

# Output

`modeldata` = the dictionary on input is augmented with the keys

  - `"geom"` = the nodal field that is the geometry
  - `"u"` = the nodal field that is the computed displacement
  - `"temp"` = the nodal field that is the temperature change   
  - `"work"` = work of the applied loads    
  - `"timing"` = dictionary with timing results    
"""
function linearstatics(modeldata::FDataDict)

    # Lists of recognized keys for the data dictionaries:
    modeldata_recognized_keys = [
        "fens",
        "regions",
        "essential_bcs",
        "traction_bcs",
        "temperature_change",
        "factorize",
    ]
    essential_bcs_recognized_keys = ["displacement", "node_list", "component"]
    traction_bcs_recognized_keys = ["femm", "traction_vector"]
    regions_recognized_keys = ["femm", "body_load"]
    temperature_change_recognized_keys = ["temperature"]

    # Extract the nodes
    fens = get(() -> error("Must get fens!"), modeldata, "fens")

    # Construct the geometry field
    geom = NodalField(fens.xyz)

    # Construct the displacement field
    u = NodalField(zeros(nnodes(geom), ndofs(geom)))

    UFT = eltype(u.values)

    # Construct the temperature field
    temp = NodalField(zeros(nnodes(geom), 1))

    modeldata["timing"] = FDataDict()

    tstart = time()
    # Apply the essential boundary conditions on the displacement field
    essential_bcs = get(modeldata, "essential_bcs", nothing)
    if (essential_bcs !== nothing)
        for j in eachindex(essential_bcs)
            ebc = essential_bcs[j]
            dcheck!(ebc, essential_bcs_recognized_keys)
            fenids = get(() -> error("Must get node list!"), ebc, "node_list")
            displacement = get(ebc, "displacement", nothing)
            u_fixed = zeros(UFT, length(fenids)) # default is  zero displacement
            if (displacement !== nothing) # if it is nonzero,
                if (typeof(displacement) <: Function) # it could be a function
                    for k in eachindex(fenids)
                        u_fixed[k] = displacement(geom.values[fenids[k], :])[1]
                    end
                else # or it could be a constant
                    fill!(u_fixed, displacement)
                end
            end
            component = get(ebc, "component", 0) # which component?
            setebc!(u, fenids[:], true, component, u_fixed)
        end
        applyebc!(u)
    end

    # Number the equations
    numberdofs!(u)           #,Renumbering_options); # NOT DONE

    modeldata["timing"]["essential_bcs"] = time() - tstart

    tstart = time()
    # Initialize the heat loads vector
    F = zeros(UFT, nalldofs(u))
    # Construct the system stiffness matrix
    K = spzeros(nalldofs(u), nalldofs(u)) # (all zeros, for the moment)
    regions = get(() -> error("Must get region list!"), modeldata, "regions")
    for i in eachindex(regions)
        region = regions[i]
        dcheck!(region, regions_recognized_keys)
        femm = region["femm"]
        # Give the  FEMM a chance  to precompute  geometry-related quantities
        femm = associategeometry!(femm, geom)
        # Add up all the conductivity matrices for all the regions
        K = K + stiffness(femm, geom, u)
    end

    tstart = time()
    # Process the traction boundary condition
    traction_bcs = get(modeldata, "traction_bcs", nothing)
    if (traction_bcs !== nothing)
        for j in eachindex(traction_bcs)
            tractionbc = traction_bcs[j]
            dcheck!(tractionbc, traction_bcs_recognized_keys)
            traction_vector = tractionbc["traction_vector"]
            if (typeof(traction_vector) <: Function)
                fi = ForceIntensity(UFT, ndofs(geom), traction_vector)
            elseif (typeof(traction_vector) <: ForceIntensity)
                fi = traction_vector
            else
                fi = ForceIntensity(traction_vector)
            end
            femm = tractionbc["femm"]
            F = F + distribloads(femm, geom, u, fi, 2)
        end
    end
    modeldata["timing"]["traction_bcs"] = time() - tstart

    tstart = time()
    # Process the thermal strain  loading
    temperature_change = get(modeldata, "temperature_change", nothing)
    if (temperature_change !== nothing)
        dcheck!(temperature_change, temperature_change_recognized_keys)
        # zero temperature change is a reasonable default
        temp = NodalField(zeros(size(fens.xyz, 1), 1))
        temperature = get(temperature_change, "temperature", nothing)
        if (temperature !== nothing) # if it is nonzero,
            if (typeof(temperature) <: Function) # it could be a function
                for k in eachindex(fens)
                    temp.values[k] = temperature(geom.values[k, :])[1]
                end
            else # or it could be a constant
                fill!(temp.values, temperature)
            end
        end
        for i in eachindex(regions)
            region = regions[i]
            femm = region["femm"]
            F = F + thermalstrainloads(femm, geom, u, temp)
        end
    end
    modeldata["timing"]["temperature_change"] = time() - tstart

    K_ff, K_fd = matrix_blocked(K, nfreedofs(u), nfreedofs(u))[(:ff, :fd)]
    F_f = vector_blocked(F, nfreedofs(u))[:f]
    U_d = gathersysvec(u, :d)

    # Loads due to the essential boundary conditions on the displacement field
    essential_bcs = get(modeldata, "essential_bcs", nothing)
    if (essential_bcs !== nothing) # there was at least one EBC applied
        F_f = F_f - K_fd * U_d
    end

    modeldata["timing"]["stiffness"] = time() - tstart

    # # Process the body load
    # body_load = get(modeldata, "body_load", nothing);
    # if (body_load  !=nothing)
    #     for j=1:length(model_data.body_load)
    #         body_load =model_data.body_load{j};
    #         femm = femm_deformation_linear (struct ('material',[],...
    #             'fes',body_load.fes,...
    #             'integration_rule',body_load.integration_rule));
    #         fi= force_intensity(struct('magn',body_load.force));
    #         F = F + distrib_loads(femm, sysvec_assembler, geom, u, fi, 3);
    #     end
    #     clear body_load fi  femm
    # end

    # # Process the nodal force boundary condition
    # if (isfield(model_data.boundary_conditions, 'nodal_force' ))
    #     for j=1:length(model_data.boundary_conditions.nodal_force)
    #         nodal_force =model_data.boundary_conditions.nodal_force{j};
    #         femm = femm_deformation_linear (struct ('material',[],...
    #             'fes',fe_set_P1(struct('conn',reshape(nodal_force.node_list,[],1))),...
    #             'integration_rule',point_rule));
    #         fi= force_intensity(struct('magn',nodal_force.force));
    #         F = F + distrib_loads(femm, sysvec_assembler, geom, u, fi, 0);
    #     end
    #     clear nodal_force fi femm
    # end

    # # Apply multi point constraints
    # if isfield(model_data,'mpc')
    #     for i=1:length(model_data.mpc)
    #         mpc =model_data.mpc{i};
    #         dofnums=0*mpc.umultipliers;# Construct an array of the degree of freedom numbers
    #         for kx=1:length(mpc.node_list)
    #             dofnums(kx)=u.dofnums(mpc.node_list(kx),mpc.dof_list(kx));
    #         end
    #         # Now call the utility function to calculate the constraint matrix
    #         [Kmpc,Fmpc]=apply_penalty_mpc(nalldofs(u),dofnums,mpc.umultipliers,0.0,mpc.penfact);
    #         K = K + Kmpc;
    #         F = F + Fmpc;
    #     end
    #     clear Kmpc Fmpc
    # end

    tstart = time()
    # Solve the system of linear algebraic equations
    U_f = K_ff \ F_f
    scattersysvec!(u, U_f)
    modeldata["timing"]["solution"] = time() - tstart

    U = gathersysvec(u, :a)

    # Update the model data
    setindex!(modeldata, geom, "geom")
    setindex!(modeldata, u, "u")
    setindex!(modeldata, temp, "temp")
    setindex!(modeldata, temp, "dT")
    setindex!(modeldata, dot(F, U) / 2, "work")
    return modeldata            # ... And return the updated model data
end

"""
    AlgoDeforLinearModule.exportdeformation(modeldata::FDataDict)

Algorithm for exporting of the deformation for visualization in Paraview.

# Argument

`modeldata` = dictionary with values for keys

  - `"fens"`  = finite element node set

  - `"regions"`  = array of region dictionaries
  - `"geom"` = geometry field
  - `"u"` = displacement field, or
  - `"us"` = array of  tuples (name, displacement field)
  - `"postprocessing"` = dictionary  with values for keys

      + `"boundary_only"` = should only the boundary of the  regions be rendered?
        Default is render the interior.
      + `"file"` = name of the  postprocessing file

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:

  - `"femm"` = finite element mmodel machine (mandatory);

# Output

`modeldata` updated with

  - `modeldata["postprocessing"]["exported"]` = array of data dictionaries, one for
    each exported file. The data is stored with the keys:

      + `"file"` - names of exported file  
      + `"field"` - nodal or elemental field
"""
function exportdeformation(modeldata::FDataDict)
    modeldata_recognized_keys = ["fens", "regions", "geom", "u", "postprocessing"]
    postprocessing_recognized_keys = ["boundary_only", "file"]
    # Defaults
    boundary_only = false
    ffile = "deformation"
    dcheck!(modeldata, modeldata_recognized_keys)

    # Let's have a look at what's been specified
    postprocessing = get(modeldata, "postprocessing", nothing)
    if (postprocessing !== nothing)
        dcheck!(postprocessing, postprocessing_recognized_keys)
        boundary_only = get(postprocessing, "boundary_only", boundary_only)
        ffile = get(postprocessing, "file", ffile)
    end

    fens = get(() -> error("Must get fens!"), modeldata, "fens")
    geom = get(() -> error("Must get geometry field!"), modeldata, "geom")
    u = get(modeldata, "u", nothing)
    UFT = eltype(u.values)
    us = get(modeldata, "us", nothing)
    if us === nothing
        us = [("u", u)]
    end

    # Export one file for each region
    modeldata["postprocessing"]["exported"] = Array{FDataDict,1}()
    regions = get(() -> error("Must get region!"), modeldata, "regions")
    for i in eachindex(regions)
        region = regions[i]
        femm = region["femm"]
        rfile = ffile * "$i" * ".vtk"
        vectors = Tuple{String,Matrix{UFT}}[]
        for ixxxx in eachindex(us)
            push!(vectors, (us[ixxxx][1], us[ixxxx][2].values))
        end
        if boundary_only
            bfes = meshboundary(femm.integdomain.fes)
            vtkexportmesh(rfile, fens, bfes; vectors = vectors)
        else
            vtkexportmesh(rfile, fens, femm.integdomain.fes; vectors = vectors)
        end
        ed = FDataDict(
            "file" => rfile,
            "field" => u,
            "region" => i,
            "type" => "displacement",
        )
        push!(modeldata["postprocessing"]["exported"], ed)
    end

    return modeldata
end

"""
    AlgoDeforLinearModule.exportstress(modeldata::FDataDict)

Algorithm for exporting of the stress for visualization in Paraview.

# Argument

`modeldata` = dictionary with values for keys

  - `"fens"`  = finite element node set

  - `"regions"`  = array of region dictionaries
  - `"geom"` = geometry field
  - `"u"` = displacement field
  - `"postprocessing"` = dictionary  with values for keys

      + `"boundary_only"` = should only the boundary of the  regions be rendered?
        Default is render the interior.
      + `"file"` = name of the  postprocessing file
      + `"quantity"` = quantity to be exported (default `:Cauchy`)
      + `"component"` = which component of the quantity?
      + `"outputcsys"` = output coordinate system
      + `"inspectormeth"` = inspector method to pass to `inspectintegpoints()`
      + `"extrap"` = method for extrapolating from the quadrature points to the nodes
        within one element

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:

  - `"femm"` = finite element mmodel machine (mandatory);

# Output

`modeldata` updated with

  - `modeldata["postprocessing"]["exported"]` = array of data dictionaries, one for
    each exported file. The data is stored with the keys:

      + `"file"` - name of exported file
      + `"field"` - nodal field
"""
function exportstress(modeldata::FDataDict)
    modeldata_recognized_keys = ["fens", "regions", "geom", "u", "dT", "postprocessing"]
    postprocessing_recognized_keys = [
        "boundary_only",
        "file",
        "quantity",
        "component",
        "outputcsys",
        "nodevalmethod",
        "reportat",
    ]
    # Defaults
    boundary_only = false
    ffile = "stress"
    dcheck!(modeldata, modeldata_recognized_keys)
    quantity = :Cauchy
    component = 1
    outputcsys = nothing
    reportat = :default
    nodevalmethod = :invdistance
    # Let's have a look at what's been specified
    postprocessing = get(modeldata, "postprocessing", nothing)
    if (postprocessing !== nothing)
        dcheck!(postprocessing, postprocessing_recognized_keys)
        boundary_only = get(postprocessing, "boundary_only", boundary_only)
        ffile = get(postprocessing, "file", ffile)
        quantity = get(postprocessing, "quantity", quantity)
        component = get(postprocessing, "component", component)
        outputcsys = get(postprocessing, "outputcsys", outputcsys)
        nodevalmethod = get(postprocessing, "nodevalmethod", nodevalmethod)
        reportat = get(postprocessing, "reportat", reportat)
    end

    fens = get(() -> error("Must get fens!"), modeldata, "fens")
    geom = get(() -> error("Must get geometry field!"), modeldata, "geom")
    u = get(() -> error("Must get displacement field!"), modeldata, "u")
    dT = get(modeldata, "dT", nothing)

    context = []
    if (outputcsys !== nothing)
        push!(context, (:outputcsys, outputcsys))
    end
    if (nodevalmethod !== nothing)
        push!(context, (:nodevalmethod, nodevalmethod))
    end
    if (reportat !== nothing)
        push!(context, (:reportat, reportat))
    end

    # Export a file for each region
    modeldata["postprocessing"]["exported"] = Array{FDataDict,1}()
    regions = get(() -> error("Must get region!"), modeldata, "regions")
    for i in eachindex(regions)
        region = regions[i]
        femm = region["femm"]
        rfile = ffile * "-" * string(quantity) * string(component) * "-region $i" * ".vtk"
        if (typeof(component) == Symbol)
            componentnum = stresscomponentmap(femm.mr)[component]
        else
            componentnum = component
        end
        componentname = length(componentnum) > 1 ? "" : "$(componentnum)"
        # Note that we are creating a field  separately for each region.  This is
        # important  for the following reason: if the regions were of different
        # materials, or if they were of the same material but with different material
        # axes orientation, averaging across the material interface  would not make
        # sense.
        if (dT !== nothing)
            fld =
                fieldfromintegpoints(femm, geom, u, dT, quantity, componentnum; context...)
        else
            fld = fieldfromintegpoints(femm, geom, u, quantity, componentnum; context...)
        end
        if boundary_only
            bfes = meshboundary(femm.integdomain.fes)
            vtkexportmesh(
                rfile,
                fens,
                bfes;
                scalars = [(string(quantity) * componentname, fld.values)],
                vectors = [("u", u.values)],
            )
        else
            vtkexportmesh(
                rfile,
                fens,
                femm.integdomain.fes;
                scalars = [(string(quantity) * componentname, fld.values)],
                vectors = [("u", u.values)],
            )
        end
        ed = FDataDict(
            "file" => rfile,
            "field" => fld,
            "region" => i,
            "type" => "nodal stress",
            "quantity" => quantity,
            "component" => component,
            "outputcsys" => outputcsys,
            "nodevalmethod" => nodevalmethod,
            "reportat" => reportat,
        )
        push!(modeldata["postprocessing"]["exported"], ed)
    end

    return modeldata
end

"""
    AlgoDeforLinearModule.exportstresselementwise(modeldata::FDataDict)

Algorithm for exporting of the elementwise stress for visualization in Paraview.

# Argument

`modeldata` = dictionary with values for keys

  - `"fens"`  = finite element node set

  - `"regions"`  = array of region dictionaries
  - `"geom"` = geometry field
  - `"u"` = displacement field
  - `"postprocessing"` = dictionary  with values for keys

      + `"boundary_only"` = should only the boundary of the  regions be rendered?
        Default is render the interior.
      + `"file"` = name of the  postprocessing file
      + `"quantity"` = quantity to be exported (default `:Cauchy`)
      + `"component"` = which component of the quantity?
      + `"outputcsys"` = output coordinate system

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:

  - `"femm"` = finite element mmodel machine (mandatory);

# Output

`modeldata` updated with

  - `modeldata["postprocessing"]["exported"]` = array of data dictionaries, one for
    each exported file. The data is stored with the keys:

      + `"file"` - name of exported file
      + `"field"` - elemental field
"""
function exportstresselementwise(modeldata::FDataDict)
    modeldata_recognized_keys = ["fens", "regions", "geom", "u", "dT", "postprocessing"]
    postprocessing_recognized_keys =
        ["boundary_only", "file", "quantity", "component", "outputcsys"]
    # Defaults
    boundary_only = false
    ffile = "stress"
    dcheck!(modeldata, modeldata_recognized_keys)
    quantity = :Cauchy
    component = 1
    outputcsys = nothing
    # Let's have a look at what's been specified
    postprocessing = get(modeldata, "postprocessing", nothing)
    if (postprocessing !== nothing)
        dcheck!(postprocessing, postprocessing_recognized_keys)
        boundary_only = get(postprocessing, "boundary_only", boundary_only)
        ffile = get(postprocessing, "file", ffile)
        quantity = get(postprocessing, "quantity", quantity)
        component = get(postprocessing, "component", component)
        outputcsys = get(postprocessing, "outputcsys", outputcsys)
    end

    fens = get(() -> error("Must get fens!"), modeldata, "fens")
    geom = get(() -> error("Must get geometry field!"), modeldata, "geom")
    u = get(() -> error("Must get displacement field!"), modeldata, "u")
    dT = get(modeldata, "dT", nothing)

    context = []
    if (outputcsys !== nothing)
        push!(context, (:outputcsys, outputcsys))
    end

    # Export a file for each region
    modeldata["postprocessing"]["exported"] = Array{FDataDict,1}()
    regions = get(() -> error("Must get region!"), modeldata, "regions")
    for i in eachindex(regions)
        region = regions[i]
        femm = region["femm"]
        rfile = ffile * "-" * string(quantity) * string(component) * "-region $i" * ".vtk"
        if (typeof(component) == Symbol)
            componentnum = stresscomponentmap(femm.mr)[component]
        else
            componentnum = component
        end
        componentname = length(componentnum) > 1 ? "" : "$(componentnum)"
        # Note that we are creating a field  separately for each region.  This is
        # important  for the following reason: if the regions were of different
        # materials, or if they were of the same material but with different material
        # axes orientation, averaging across the material interface  would not make
        # sense.
        if (dT !== nothing)
            fld = elemfieldfromintegpoints(
                femm,
                geom,
                u,
                dT,
                quantity,
                componentnum;
                context...,
            )
        else
            fld =
                elemfieldfromintegpoints(femm, geom, u, quantity, componentnum; context...)
        end
        if boundary_only
            bfes = meshboundary(femm.integdomain.fes)
            vtkexportmesh(
                rfile,
                fens,
                bfes;
                scalars = [(string(quantity) * componentname, fld.values)],
                vectors = [("u", u.values)],
            )
        else
            vtkexportmesh(
                rfile,
                fens,
                femm.integdomain.fes;
                scalars = [(string(quantity) * componentname, fld.values)],
                vectors = [("u", u.values)],
            )
        end
        ed = FDataDict(
            "file" => rfile,
            "field" => fld,
            "region" => i,
            "type" => "elemental stress",
            "quantity" => quantity,
            "component" => component,
            "outputcsys" => outputcsys,
        )
        push!(modeldata["postprocessing"]["exported"], ed)
    end

    return modeldata
end

"""
    AlgoDeforLinearModule.modal(modeldata::FDataDict)

Modal (free-vibration) analysis solver.

# Argument

`modeldata` = dictionary with values for keys

  - `"fens"`  = finite element node set
  - `"regions"`  = array of region dictionaries
  - `"essential_bcs"` = array of essential boundary condition dictionaries

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:

  - `"femm"` = finite element mmodel machine (mandatory);

For essential boundary conditions (optional) each dictionary
would hold

  - `"displacement"` = fixed (prescribed) displacement (scalar): only zero
    displacement is  allowed for modal analysis.
  - `"component"` = which component is prescribed  (1, 2, 3)?
  - `"node_list"` = list of nodes on the boundary to which the condition applies
    (mandatory)

Control parameters:

  - `"neigvs"` = number of eigenvalues/eigenvectors to compute
  - `"omega_shift"`= angular frequency shift for mass shifting
  - `"use_lumped_mass"` = true or false?  (Default is false: consistent mass)

# Output

`modeldata`= the dictionary on input is augmented with

  - `"geom"` = the nodal field that is the geometry
  - `"u"` = the nodal field that is the computed displacement
  - `"neigvs"` = Number of computed eigenvectors
  - `"W"` = Computed eigenvectors, neigvs columns
  - `"omega"` =  Computed angular frequencies, array of length neigvs    # For multi point constraints (MPC) (optional):
  - `"raw_eigenvalues"` = Raw computed eigenvalues    # model_data.mpc= cell array of structs, each for one MPC.
"""
function modal(modeldata::FDataDict)

    # For multi point constraints (MPC) (optional):
    # model_data.mpc= cell array of structs, each for one MPC.
    #      mpc.node_list = list of node numbers involved in the MPC,
    #      mpc.dof_list= numbers of degrees of freedom for the nodes above,
    #      mpc.umultipliers=multipliers for the nodes above,
    #      mpc.penfact=the penalty factor to multiply  the constraint matrix,
    #          The MPC looks like this: sum_i m_i u_{dof(i),node(i)} =0
    #          where m_i is the multiplier.

    # Lists of recognized keys for the data dictionaries:
    modeldata_recognized_keys =
        ["fens", "regions", "essential_bcs", "neigvs", "omega_shift", "use_lumped_mass"]
    essential_bcs_recognized_keys = ["displacement", "node_list", "component"]
    regions_recognized_keys = ["femm", "femm_stiffness", "femm_mass", "body_load"]

    neigvs = get(modeldata, "neigvs", 7) # Number of eigenvalues

    omega_shift = get(modeldata, "omega_shift", 0.0) # Mass shifting

    use_factorization = get(modeldata, "use_factorization", false) # Factorization?

    use_lumped_mass = get(modeldata, "use_lumped_mass", false) # Lumped mass?

    # Extract the nodes
    fens = get(() -> error("Must get fens!"), modeldata, "fens")

    # Construct the geometry field
    geom = NodalField(fens.xyz)

    # Construct the displacement field
    u = NodalField(zeros(nnodes(geom), ndofs(geom)))
    UFT = eltype(u.values)

    # Apply the essential boundary conditions on the displacement field
    essential_bcs = get(modeldata, "essential_bcs", nothing)
    if (essential_bcs !== nothing)
        for j in eachindex(essential_bcs)
            ebc = essential_bcs[j]
            dcheck!(ebc, essential_bcs_recognized_keys)
            fenids = get(() -> error("Must get node list!"), ebc, "node_list")
            displacement = get(ebc, "displacement", nothing)
            u_fixed = zeros(UFT, length(fenids)) # only zero displacement accepted
            component = get(ebc, "component", 0) # which component?
            setebc!(u, fenids[:], true, component, u_fixed)
        end
        applyebc!(u)
    end

    # Number the equations
    numberdofs!(u)           #,Renumbering_options); # NOT DONE

    # Construct the system stiffness matrix
    K = spzeros(nalldofs(u), nalldofs(u)) # (all zeros, for the moment)
    regions = get(() -> error("Must get region list!"), modeldata, "regions")
    for i in eachindex(regions)
        region = regions[i]
        dcheck!(region, regions_recognized_keys)
        if "femm_stiffness" in keys(region)
            femm = region["femm_stiffness"]
        else
            femm = get(() -> error("Must get femm or femm_stiffness!"), region, "femm")
        end
        # Give the  FEMM a chance  to precompute  geometry-related quantities
        femm = associategeometry!(femm, geom)
        # Add up all the stiffness matrices for all the regions
        K = K + stiffness(femm, geom, u)
    end

    # Construct the system mass matrix
    M = spzeros(nalldofs(u), nalldofs(u)) # (all zeros, for the moment)
    regions = get(() -> error("Must get region list!"), modeldata, "regions")
    for i in eachindex(regions)
        region = regions[i]
        dcheck!(region, regions_recognized_keys)
        if "femm_mass" in keys(region)
            femm = region["femm_mass"]
        else
            femm = get(() -> error("Must get femm or femm_mass!"), region, "femm")
        end
        # Give the  FEMM a chance  to precompute  geometry-related quantities
        femm = associategeometry!(femm, geom)
        # Add up all the mass matrices for all the regions
        M = M + mass(femm, geom, u)
    end

    # Options for the eigenproblem solution

    # Solve
    # if (~ use_factorization )
    #     # This is one way of solving the eigenvalue problem, just pass the matrices
    #     [W,Omega]= eigs(K+omega_shift*M, M, neigvs, 'SM', evopts);
    # else
    # This form uses the factorized matrix and has the potential of being much faster
    # Factorize the left-hand side matrix for efficiency (Choleski)
    # [mA,status] = chol(K+omega_shift*M,'lower');#,'vector',prm
    # if ( status ~= 0 ) error('Choleski factorization failed'), end
    # clear K; # Not needed anymore
    # mAt= mA';
    # [W,Omega]= eigs(@(bv)mAt\(mA\bv), nalldofs(u), M, neigvs, 'SM', evopts);
    #          [W,Omega]= eigen(full(K+omega_shift*M), full(M));

    K_ff = matrix_blocked(K, nfreedofs(u), nfreedofs(u))[:ff]
    M_ff = matrix_blocked(M, nfreedofs(u), nfreedofs(u))[:ff]

    d, v, nconv = eigs(Symmetric(K_ff + omega_shift * M_ff), Symmetric(M_ff); nev = neigvs, which = :SM, explicittransform=:none)
    #    Subtract the mass-shifting Angular frequency
    broadcast!(+, d, d, -omega_shift)

    modeldata["raw_eigenvalues"] = d
    # Better make sure the eigenvalues make sense: they should be all real.
    if any(imag(d) .!= 0.0)
        d = real.(d)
    end
    if any(real(d) .< 0.0)
        d = abs.(d)
    end
    d = real.(d)
    #    Sort  the angular frequencies by magnitude.  Make sure all
    #    imaginary parts of the eigenvalues are removed.
    ix = sortperm(d)

    # Update the model data: store geometry
    modeldata["geom"] = geom
    # Store the displacement field
    modeldata["u"] = u
    # Number of computed eigenvectors
    modeldata["neigvs"] = length(d)
    #  Computed eigenvectors: we are ignoring the imaginary part here
    #  because the modal analysis is presumed to have been performed for
    #  an undamped structure
    modeldata["W"] = real(v[:, ix])
    #  Computed angular frequencies
    modeldata["omega"] = sqrt.(d[ix])
    return modeldata
end

"""
    AlgoDeforLinearModule.exportmode(modeldata::FDataDict)

Algorithm for exporting of the mmode shape for visualization in Paraview.

# Argument

`modeldata` = dictionary with values for keys

  - `"fens"`  = finite element node set

  - `"regions"`  = array of region dictionaries
  - `"geom"` = geometry field
  - `"u"` = displacement field
  - `"W"` = Computed free-vibration eigenvectors, `neigvs` columns
  - `"omega"` =  Computed free-vibration angular frequencies, array of length `neigvs`
  - `"postprocessing"` = dictionary  with values for keys

      + `"boundary_only"` = should only the boundary of the  regions be rendered?
        Default is render the interior.
      + `"file"` = name of the  postprocessing file
      + `"mode"` = which mode should be visualized?
      + `"component"` = which component of the quantity?
      + `"outputcsys"` = output coordinate system

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:

  - `"femm"` = finite element mmodel machine (mandatory);

# Output

`modeldata` updated with

  - `modeldata["postprocessing"]["exported"]` = see `exportdeformation()`
"""
function exportmode(modeldata::FDataDict)
    modeldata_recognized_keys =
        ["fens", "regions", "geom", "u", "omega", "W", "postprocessing"]
    postprocessing_recognized_keys = ["boundary_only", "file", "mode"]
    mode = 1
    dcheck!(modeldata, modeldata_recognized_keys)

    # Let's have a look at what's been specified
    postprocessing = get(modeldata, "postprocessing", nothing)
    if (postprocessing !== nothing)
        dcheck!(postprocessing, postprocessing_recognized_keys)
        mode = get(postprocessing, "mode", mode)
    end

    omega = modeldata["omega"]

    # Scatter the desired mode
    W = modeldata["W"]
    if typeof(mode) <: Int
        @assert 0 < mode <= length(omega) "Invalid mode number $mode"
        scattersysvec!(modeldata["u"], W[:, mode])
    else
        us = Tuple{String,AbstractField}[]
        u = modeldata["u"]
        for ixxxx in mode
            @assert 0 < ixxxx <= length(omega) "Invalid mode number $ixxxx"
            scattersysvec!(u, W[:, ixxxx])
            push!(us, ("mode_$(ixxxx)", deepcopy(u)))
        end
        modeldata["us"] = us
    end

    return exportdeformation(modeldata)
end

end
