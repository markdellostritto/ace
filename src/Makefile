##################################################################################
# EXTERNAL LIBRARY PATHS
##################################################################################

EIGEN = /usr/local/include/eigen-3.4-rc1/ # eigen library
SRC  := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
INC_LIST = $(EIGEN) $(SRC) 
INC = $(foreach d, $(INC_LIST), -I$d)

##################################################################################
# GENERAL SETTINGS
##################################################################################
COMP=gnu

# Include global settings.
include Makefile.$(COMP)

##################################################################################
# OBJECT FILES
##################################################################################

util = compiler.o time.o
struc = structure.o cell.o state.o atom.o neighbor.o 
math = special.o reduce.o random.o eigen.o corr.o
chem = units.o ptable.o 
str = string.o print.o token.o
fstruc = file_struc.o vasp_struc.o xyz_struc.o cube_struc.o qe_struc.o
mem = serialize.o map.o
calc = calc.o calc_factory.o calc_lj_cut.o calc_coul_cut.o calc_coul_long.o calc_cgem_cut.o calc_cgem_long.o
opt = loss.o stop.o algo.o objective.o decay.o iter.o
constraint = constraint.o constraint_factory.o constraint_freeze.o 
sim = engine.o integrator.o job.o
kspaces = kspace.o kspace_coul.o
thread = comm.o dist.o
neural = batch.o nn.o 
basis = basis.o basis_radial.o basis_angular.o
nnp = nnh.o nnp.o cutoff.o 

##################################################################################
# OBJECT FILE MAKE RULES
##################################################################################

# ==== struc ====

structure.o: struc/structure.cpp 
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c struc/structure.cpp
cell.o: struc/cell.cpp 
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c struc/cell.cpp
atom.o: struc/atom.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c struc/atom.cpp
state.o: struc/state.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c struc/state.cpp
neighbor.o: struc/neighbor.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c struc/neighbor.cpp
grid.o: struc/grid.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c struc/grid.cpp

# ==== math ====

special.o: math/special.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c math/special.cpp
reduce.o: math/reduce.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c math/reduce.cpp
random.o: math/random.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c math/random.cpp
eigen.o: math/eigen.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c math/eigen.cpp
corr.o: math/corr.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c math/corr.cpp

# ==== chem ====

units.o: chem/units.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c chem/units.cpp
ptable.o: chem/ptable.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c chem/ptable.cpp
alias.o: chem/alias.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c chem/alias.cpp

# ==== string ====

string.o: str/string.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c str/string.cpp
print.o: str/print.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c str/print.cpp
token.o: str/token.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c str/token.cpp

# ==== utility ====

compiler.o: util/compiler.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c util/compiler.cpp
time.o: util/time.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c util/time.cpp

# ==== memory ====

serialize.o: mem/serialize.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c mem/serialize.cpp
map.o: mem/map.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c mem/map.cpp

# ==== thread ====

dist.o: thread/dist.cpp
	$(CXX_THREAD) $(CXX_FLAGS) $(INC) -c thread/dist.cpp
comm.o: thread/comm.cpp
	$(CXX_THREAD) $(CXX_FLAGS) $(INC) -c thread/comm.cpp

# ==== molecular dynamics ====

job.o: sim/job.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/job.cpp
integrator.o: sim/integrator.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/integrator.cpp
engine.o: sim/engine.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/engine.cpp

# ==== format ====

format.o: format/format.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c format/format.cpp
file_struc.o: format/file_struc.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c format/file_struc.cpp
vasp_struc.o: format/vasp_struc.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c format/vasp_struc.cpp
xyz_struc.o: format/xyz_struc.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c format/xyz_struc.cpp
cube_struc.o: format/cube_struc.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c format/cube_struc.cpp
qe_struc.o: format/qe_struc.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c format/qe_struc.cpp

# ==== machine learning ====

nn.o: ml/nn.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c ml/nn.cpp
batch.o: ml/batch.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c ml/batch.cpp
nn_train.o: ml/nn_train.cpp
	$(CXX_THREAD) $(CXX_FLAGS) $(INC) -c ml/nn_train.cpp
pca.o: ml/pca.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c ml/pca.cpp

# ==== basis ====

basis.o: nnp/basis.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c nnp/basis.cpp
basis_radial.o: nnp/basis_radial.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c nnp/basis_radial.cpp
basis_angular.o: nnp/basis_angular.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c nnp/basis_angular.cpp

# ==== nnp ====

nnh.o: nnp/nnh.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c nnp/nnh.cpp
nnp.o: nnp/nnp.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c nnp/nnp.cpp
cutoff.o: nnp/cutoff.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c nnp/cutoff.cpp

# ==== optimization ====

loss.o: opt/loss.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c opt/loss.cpp
stop.o: opt/stop.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c opt/stop.cpp
algo.o: opt/algo.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c opt/algo.cpp
decay.o: opt/decay.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c opt/decay.cpp
objective.o: opt/objective.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c opt/objective.cpp
iter.o: opt/iter.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c opt/iter.cpp

# ==== calc ====

calc.o: sim/calc.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/calc.cpp
calc_factory.o: sim/calc_factory.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/calc_factory.cpp
calc_lj_cut.o: sim/calc_lj_cut.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/calc_lj_cut.cpp
calc_coul_cut.o: sim/calc_coul_cut.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/calc_coul_cut.cpp
calc_coul_long.o: sim/calc_coul_long.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/calc_coul_long.cpp
calc_cgem_cut.o: sim/calc_cgem_cut.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/calc_cgem_cut.cpp
calc_cgem_long.o: sim/calc_cgem_long.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/calc_cgem_long.cpp

# ==== constraint ====

constraint.o: sim/constraint.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/constraint.cpp
constraint_factory.o: sim/constraint_factory.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/constraint_factory.cpp
constraint_freeze.o: sim/constraint_freeze.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/constraint_freeze.cpp

# ==== kspaces ====

kspace.o: sim/kspace.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/kspace.cpp
kspace_coul.o: sim/kspace_coul.cpp
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -c sim/kspace_coul.cpp

##################################################################################
# TARGETS
##################################################################################

clean: 
	rm *.o

# ==== UNIT TESTS ====

#** math **
test_math_special: special.o print.o reduce.o
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -o test_math_special.exe math/test_math_special.cpp special.o print.o reduce.o
test_math_reduce: reduce.o print.o
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -o test_math_reduce.exe math/test_math_reduce.cpp reduce.o print.o
# chem
test_units: units.o print.o
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -o test_units.exe chem/test_units.cpp units.o print.o
# string
test_token: $(str)
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -o test_token.exe str/test_token.cpp $(str)

# ==== MOLECULAR DYNAMICS ====
md: $(str) $(math) $(chem) $(struc) $(fstruc) $(mem) $(kspaces) $(calc) $(constraint) $(sim) $(util) format.o grid.o
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -o md.exe sim/md.cpp $(str) $(math) $(chem) $(struc) $(fstruc) $(mem) $(kspaces) $(calc) $(constraint) $(sim) $(util) format.o grid.o 
fit_cgem: $(str) $(math) $(chem) $(struc) $(fstruc) $(mem) $(calc) $(kspaces) $(constraint) $(sim) format.o grid.o -lnlopt
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -o fit_cgem.exe sim/fit_cgem.cpp $(str) $(math) $(chem) $(struc) $(fstruc) $(mem) $(calc) $(kspaces) $(constraint) $(sim) format.o grid.o -lnlopt

# ==== NEURAL NETWORK ====
nn_fit: $(neural) $(str) $(math) $(util) $(thread) $(opt) serialize.o compiler.o nn_train.o
	$(CXX_THREAD) $(CXX_FLAGS) $(INC) -o nn_fit.exe ml/nn_fit.cpp $(neural) $(str) $(math) $(util) $(thread) $(opt) serialize.o compiler.o nn_train.o

# ==== NNP ====
test_nnp_cutoff: $(mem) $(str) $(math) basis.o cutoff.o
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -o test_nnp_cutoff.exe nnp/test_nnp_cutoff.cpp $(mem) $(str) $(math) basis.o cutoff.o
test_basis_radial: $(mem) $(str) $(math) basis.o basis_radial.o cutoff.o
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -o test_basis_radial.exe nnp/test_basis_radial.cpp $(mem) $(str) $(math) basis.o basis_radial.o cutoff.o
test_basis_angular: $(mem) $(str) $(math) basis.o basis_angular.o cutoff.o
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -o test_basis_angular.exe nnp/test_basis_angular.cpp $(mem) $(str) $(math) basis.o basis_angular.o cutoff.o
test_nnp_symm: $(basis) $(mem) $(struc) $(str) $(math) nnh.o nnp.o nn.o cutoff.o 
	$(CXX_SERIAL) $(CXX_FLAGS) $(INC) -o test_nnp_symm.exe nnp/test_nnp_symm.cpp $(basis) $(mem) $(struc) $(str) $(math) nnh.o nnp.o nn.o cutoff.o 
test_nnp_force: $(mem) $(struc) $(fstruc) $(str) $(math) $(chem) $(basis) $(nnp) $(neural) $(opt) $(thread) grid.o
	$(CXX_THREAD) $(CXX_FLAGS) $(INC) -o test_nnp_force.exe nnp/test_nnp_force.cpp $(mem) $(struc) $(fstruc) $(str) $(math) $(chem) $(basis) $(nnp) $(neural) $(opt) $(thread) grid.o
nnptefr2: $(nnp) $(str) $(math) $(mem) $(thread) $(basis) $(struc) $(fstruc) $(chem) $(kspaces) $(opt) $(neural) $(util) $(calc) $(constraint) format.o alias.o pca.o grid.o engine.o integrator.o 
	$(CXX_THREAD) $(CXX_FLAGS) $(INC) -o nnptefr2.exe nnp/nnptefr2.cpp $(nnp) $(str) $(math) $(mem) $(thread) $(basis) $(struc) $(fstruc) $(chem) $(kspaces) $(opt) $(neural) $(util) $(calc) $(constraint) format.o alias.o pca.o grid.o engine.o integrator.o 
nnptefr: $(nnp) $(str) $(math) $(mem) $(thread) $(basis) $(struc) $(fstruc) $(chem) $(kspaces) $(opt) $(neural) $(util) $(calc) $(constraint) format.o alias.o pca.o grid.o engine.o integrator.o
	$(CXX_THREAD) $(CXX_FLAGS) $(INC) -o nnptefr.exe nnp/nnptefr.cpp $(nnp) $(str) $(math) $(mem) $(thread) $(basis) $(struc) $(fstruc) $(chem) $(kspaces) $(opt) $(neural) $(util) $(calc) $(constraint) format.o alias.o pca.o grid.o engine.o integrator.o 

