FC      = gfortran
FFLAGS  = -O2 -fimplicit-none -Wall -Wextra -Jbuild -Ibuild -g -fcheck=bounds
SRC_DIR = src
BUILD   = build
BIN     = vertex3d

# Compile order matters for module dependencies.
SOURCES = \
  $(SRC_DIR)/mod_kinds.f90 \
  $(SRC_DIR)/mod_parameters.f90 \
  $(SRC_DIR)/mod_data.f90 \
  $(SRC_DIR)/mod_geometry.f90 \
  $(SRC_DIR)/mod_topology.f90 \
  $(SRC_DIR)/mesh_init.f90 \
  $(SRC_DIR)/Force.f90 \
  $(SRC_DIR)/Langevin_update.f90 \
  $(SRC_DIR)/T1_transition.f90 \
  $(SRC_DIR)/T2_transition.f90 \
  $(SRC_DIR)/T4_transition.f90 \
  $(SRC_DIR)/Cell_Division.f90 \
  $(SRC_DIR)/io_module.f90 \
  $(SRC_DIR)/sanity_checks.f90 \
  $(SRC_DIR)/main.f90

OBJECTS = $(patsubst $(SRC_DIR)/%.f90,$(BUILD)/%.o,$(SOURCES))

all: $(BIN)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: $(SRC_DIR)/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BIN): $(OBJECTS)
	$(FC) $(FFLAGS) -o $(BIN) $(OBJECTS)

run: $(BIN)
	./$(BIN)

clean:
	rm -rf $(BUILD) $(BIN)

distclean: clean
	rm -rf data

.PHONY: all run clean distclean
