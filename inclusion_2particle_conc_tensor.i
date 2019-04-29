# modified modelB_2particle including bodyforce irrad from irrad_Hex_poly.i deck,
# model_C_example conserved variables, and TensorMechanics from inclusion.i deck
Body_Force = 10

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 64
  ny = 64
  xmin = 0
  xmax = 256
  ymin = 0
  ymax = 256
  uniform_refine = 2
  elem_type = QUAD4
[]

[GlobalParams]
  block = 0
  op_num = 3
  var_name_base = eta
[]

[Variables]
  [./c]
  [../]
  [./w]
  [../]
  # from irrad_Hex_poly.i
  [./eta0]
  [../]
  [./eta1]
  [../]
  [./eta2]
  [../]
  [./etab]
  [../]
  ##########################################################
  # FROM ComputeConcentrationDependentElasticityTensor deck
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
##############################################################
[]
# aux varaibles to track the free energy change (must decrease with time)
[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  # FROM ComputeConcentrationDependentElasticityTensor deck
  [./s11_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  #################################################################
  # from inclusion.i deck
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
  [./fel_an]
    order = CONSTANT
    family = MONOMIAL
  [../]
#############################################################
  # the chemical potential gradients
  [./dwdx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dwdy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # the flux
  [./jx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./jy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./j_tot]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[ICs]
  # from irrad_Hex_poly.i
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
  [./IC_c]
    type = SpecifiedSmoothCircleIC
    variable = c
    x_positions = '50 83'
    y_positions = '64 64'
    z_positions = '0 0'
    radii = '10 15'
    invalue = 1.0
    outvalue = 0.00786 # equilibrium solubility with particle of radius 15
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
### newly added for tensor stuff
  [./left]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./bottom]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  ## missfit strain workshop 9
  ## iradiation sent by email
#####################################
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
    value = ${Body_Force}  # dose rate
  [../]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
[]

[AuxKernels]
  [./bnds]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  # FROM ComputeConcentrationDependentElasticityTensor deck
  [./matl_s11]
     type = RankTwoAux
     variable = s11_aux
     rank_two_tensor = stress
     index_i = 0
     index_j = 0
   [../]
   # from inclusion.i
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
   [./matl_e11_an]
     type = RankTwoAux
     rank_two_tensor = strain_an
     index_i = 0
     index_j = 0
     variable = e11_an
   [../]
   [./matl_e12_an]
     type = RankTwoAux
     rank_two_tensor = strain_an
     index_i = 0
     index_j = 1
     variable = e12_an
   [../]
   [./matl_e22_an]
     type = RankTwoAux
     rank_two_tensor = strain_an
     index_i = 1
     index_j = 1
     variable = e22_an
   [../]
   [./matl_fel_an]
     type = MaterialRealAux
     variable = fel_an
     property = fel_an_mat
   [../]
###############################################################
   [./dwdx]
     type = VariableGradientComponent
     variable = dwdx
     gradient_variable = w
     component = x
   [../]
   [./dwdy]
     type = VariableGradientComponent
     variable = dwdy
     gradient_variable = w
     component = y
   [../]
   [./jx]
     type = ParsedAux
     variable = jx
     args = 'dwdx'
     function = '-1.0*dwdx'
   [../]
   [./jy]
     type = ParsedAux
     variable = jy
     args = 'dwdy'
     function = '-1.0*dwdy'
   [../]
   [./j_tot]
     type = ParsedAux
     variable = j_tot
     args = 'jx jy'
     function = 'sqrt(jx^2+jy^2)'
   [../]
[]

[UserObjects]
  [./hex_ic]
    type = PolycrystalHex
    # dont change to jp otherwise it will crash
# (GlobalParams/op_num):
# Unable to find a valid grain to op coloring, Make sure you have created enough variables to hold a
# valid polycrystal initial condition (no grains represented by the same variable should be allowed to
# touch, ~8 for 2D, ~25 for 3D)?

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

[Materials]
  [./FreeEng]
    type = DerivativeParsedMaterial
    args = 'c eta1 eta0 etab eta2'
    constant_names = 'a s g b'
    constant_expressions = '20.0 1.50 1.50 2.0'
    function ='sumeta:=eta0^2+eta1^2+eta2^2+etab^2;f0:=etab^2/sumeta;f2:=a*(c-f0)^2;f3:=1.0/4.0+1.0/4.0*eta1^4-1.0/2.0*eta1^2+1.0/4.0*eta0^4-1.0/2.0*eta0^2+1.0/4.0*eta2^4-1.0/2.0*eta2^2+1.0/4.0*etab^4-1.0/2.0*etab^2;f4:=s*(etab^2*eta1^2+etab^2*eta0^2+etab^2*eta2^2)+g*(eta1^2*eta0^2+eta1^2*eta2^2+eta0^2*eta2^2);f:=f2+b*(f3+f4);f'
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
  [./var_dependence]
    type = DerivativeParsedMaterial
    block = 0
    function = 0.005*c^2
    args = c
    outputs = exodus
    output_properties = 'var_dep'
    f_name = var_dep
    enable_jit = true
    derivative_order = 2
  [../]
  [./eigenstrain]
    type = ComputeVariableEigenstrain
    block = 0
    eigen_base = '1 1 0 0 0 0'
    prefactor = var_dep
    args = c
    eigenstrain_name = eigenstrain
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = 0
    displacements = 'disp_x disp_y'
    eigenstrain_names = eigenstrain
  [../]
  [./analytical]
    type = InclusionProperties
    a = 0.1
    b = 0.1
    lambda = 1
    mu = 1
    misfit_strains = '0.005 0.005'
    strain_name = strain_an
    stress_name = stress_an
    energy_name = fel_an_mat
  [../]
[]

[Postprocessors]
  [./_dt]
    # time step
    type = TimestepSize
  [../]
  [./nonlinear_its]
    # number of nonlinear iterations at each timestep
    type = NumNonlinearIterations
  [../]
  [./ElementInt_c]
    type = ElementIntegralVariablePostprocessor
    variable = c
  [../]
  [./precipitate_area]      # Area of precipitate
    type = ElementIntegralMaterialProperty
    mat_prop = prec_indic
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

[VectorPostprocessors]
  # The numerical values of the variables/auxvariables across the centerline
  [./line_values]
   type =  LineValueSampler
    start_point = '0 64 0'
    end_point = '256 64 0'
    variable = 'c w j_tot'
    num_points = 257
    sort_by =  id
    execute_on = 'TIMESTEP_END'
  [../]
[]

[Preconditioning]
  [./SMP] # to produce the complete perfect Jacobian
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
  l_max_its = 35
  l_tol = 1.0e-3
  nl_rel_tol = 1.0e-7
  start_time = 0.0
  num_steps = 200
  nl_abs_tol = 1e-8
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    growth_factor = 1.5
    cutback_factor = 0.5
    #optimal_iterations = 6
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
  #interval = 1
  file_base = 'inclusion_2_particle_src_${Body_Force}'
[]
