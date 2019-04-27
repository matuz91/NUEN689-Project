# This simulates the nucleation and growth of voids under irradiation in a ploycrystalline material
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 64
  ny = 64
  nz = 0
  xmin = 0
  xmax = 128
  ymin = 0
  ymax = 128
  zmax = 0
  uniform_refine = 1
  elem_type = QUAD4
[]

[GlobalParams]
  block = 0
  op_num = 3
  var_name_base = eta
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./c]
  [../]
  [./w]
  [../]
  [./eta0]
  [../]
  [./eta1]
  [../]
  [./eta2]
  [../]
  [./etab]
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./s11_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s12_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s22_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s11_an]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s12_an]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s22_an]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e11_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e12_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e22_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e11_an]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e12_an]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e22_an]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[UserObjects]
  [./hex_ic]
    type = PolycrystalHex
    coloring_algorithm = bt
    grain_num = 36
    x_offset = 0.0
    output_adjacency_matrix = false
  [../]
  [./normal_noise]
    type = ConservedNormalNoise
  [../]
  [./grain_tracker]
    type = GrainTracker
    threshold = 0.2
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_begin'
    flood_entity_type = ELEMENTAL
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = hex_ic
    [../]
  [../]
  [./etab] #variable representing the bubble/pore (or second-phase particle)
  type = RandomIC
  variable = etab
   max = 0.003
   min = 0.001
  [../]

  [./IC_c] #vacancy/solute concentration
  type = RandomIC
  variable = c
   max = 3e-6
   min = 1e-6
  [../]
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y'
    [../]
  [../]
  [./bottom_y]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./left_x]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
[]

[Kernels]
  [./c_dot]
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
  [./c_res]
    type = SplitCHParsed
    variable = c
    f_name = F
    kappa_name = kappa_c
    w = w
    args = 'etab eta1 eta0 eta2'
  [../]
  [./w_res]
    type = SplitCHWRes
    variable = w
    mob_name = M
  [../]
  [./AC0_bulk]
    type = AllenCahn
    variable = eta0
    f_name = F
    args = 'c eta1 etab eta2'
  [../]
  [./AC0_int]
    type = ACInterface
    variable = eta0

  [../]
  [./e0_dot]
    type = TimeDerivative
    variable = eta0
  [../]
  [./AC1_bulk]
    type = AllenCahn
    variable = eta1
    f_name = F
    args = 'c etab eta0 eta2'
  [../]
  [./AC1_int]
    type = ACInterface
    variable = eta1
  [../]
  [./e1_dot]
    type = TimeDerivative
    variable = eta1
  [../]
  [./AC2_bulk]
    type = AllenCahn
    variable = eta2
    f_name = F
    args = 'c etab eta0 eta1'
  [../]
  [./AC2_int]
    type = ACInterface
    variable = eta2
  [../]
  [./e2_dot]
    type = TimeDerivative
    variable = eta2
  [../]
  [./ACb_bulk]
    type = AllenCahn
    variable = etab
    f_name = F
    args = 'c eta1 eta0 eta2'
  [../]
  [./ACb_int]
    type = ACInterface
    variable = etab
  [../]
  [./eb_dot]
    type = TimeDerivative
    variable = etab
  [../]
  [./conserved_langevin] # Adds fluctuations in the flux without changing the total concentration
    type = ConservedLangevinNoise
    amplitude = 0.01
    variable = w
    noise = normal_noise
  [../]
  [./irrad_source] # Source term for irradiation
    type = BodyForce
    variable = w
    value = 1e-3  # dose rate
  [../]
  [./TensorMechanics]
  [../]
[]

[AuxKernels]
  [./bnds]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./matl_s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = s11_aux
  [../]
  [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = s12_aux
  [../]
  [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = s22_aux
  [../]
  [./matl_e11]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 0
    variable = e11_aux
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 1
    variable = e12_aux
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 1
    variable = e22_aux
  [../]
[]

[Materials]
  [./FreeEng]
    type = DerivativeParsedMaterial
    args = 'c eta1 eta0 etab eta2'
    block = 0
    constant_names = 'a s g b'
    constant_expressions = '20.0 1.50 1.50 2.0'
    function ='sumeta:=eta0^2+eta1^2+eta2^2+etab^2;f0:=etab^2/sumeta;f2:=a*(c-f0)^2;f3:=1.0/4.0+1.0/4.0*eta1^4-1.0/2.0*eta1^2+1.0/4.0*eta0^4-1.0/2.0*eta0^2+1.0/4.0*eta2^4-1.0/2.0*eta2^2+1.0/4.0*etab^4-1.0/2.0*etab^2;f4:=s*(etab^2*eta1^2+etab^2*eta0^2+etab^2*eta2^2)+g*(eta1^2*eta0^2+eta1^2*eta2^2+eta0^2*eta2^2);f:=f2+b*(f3+f4);f'
    #outputs = exodus
    #output_properties = 'var_dep'
    #f_name = var_dep
    enable_jit = true
    derivative_order = 2
  [../]
  [./const]
    type = GenericConstantMaterial
    prop_names = 'kappa_c kappa_op L M'
    prop_values = '0.0 0.50 1.0 1.0'
  [../]
  [./precipitate_indicator]  # Returns 1 if inside the precipitate
    type = ParsedMaterial
    f_name = prec_indic
    args = c
    function = if(c>0.9,1.0,0)
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    C_ijkl = '1 1'
    fill_method = symmetric_isotropic
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = 0
  [../]
  [./eigenstrain]
    type = ComputeVariableEigenstrain
    block = 0
    eigen_base = '1 1 0 0 0 0'
    #prefactor = var_dep
    args = c
    eigenstrain_name = eigenstrain
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = 0
    displacements = 'disp_x disp_y'
    eigenstrain_names = eigenstrain
  [../]
[]

[Postprocessors]
  [./precipitate_area]      # Area of precipitate
    type = ElementIntegralMaterialProperty
    mat_prop = prec_indic
    execute_on = 'TIMESTEP_END'
  [../]
  [./_dt]
    # time step
    type = TimestepSize
  [../]
  [./ElementInt_c]
    type = ElementIntegralVariablePostprocessor
    variable = c
  [../]
  [./num_precipitates]          # Number of precipitates
    type = FeatureFloodCount
    variable = c
    threshold = 0.9
  [../]
  [./stress_s11_aux]
    type = ElementAverageValue
    variable = s11_aux
  [../]
  [./stress_s12_aux]
    type = ElementAverageValue
    variable = s12_aux
  [../]
  [./stress_s22_aux]
    type = ElementAverageValue
    variable = s22_aux
  [../]
  [./strain_e11_aux]
    type = ElementAverageValue
    variable = e11_aux
  [../]
  [./strain_e12_aux]
    type = ElementAverageValue
    variable = e12_aux
  [../]
  [./strain_e22_aux]
    type = ElementAverageValue
    variable = e22_aux
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  nl_max_its = 15
  scheme = bdf2
  solve_type = NEWTON
  petsc_options_iname = -pc_type
  petsc_options_value = asm
  l_max_its = 30
  l_tol = 1.0e-3
  nl_rel_tol = 1.0e-8
  start_time = 0.0
  num_steps = 200
  nl_abs_tol = 1e-8
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 6
  [../]
  [./Adaptivity]
    refine_fraction = 0.5
    coarsen_fraction = 0.01
    max_h_level = 3
    initial_adaptivity = 2
  [../]
[]

[Outputs]
    exodus = true
    csv = true
    interval = 1
[]
