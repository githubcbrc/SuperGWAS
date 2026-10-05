CXX      ?= g++
CXXFLAGS  = -O3 -march=x86-64 -mtune=generic -std=c++17
EIGEN     = -I/usr/include/eigen3
VENDOR    = -I$(SRCDIR)

SRCDIR    = src
BINDIR    = bin

# Final binary list (14 total)
BINS = sample_kmers filter_kmers compute_structure convert_matrix kgwas \
       build_index map_to_ref map_merge_filter \
       compute_structure_all kgwas_2d map_2d_pairs \
       compute_structure_2w compute_projections_2w kgwas_2w

TARGETS = $(addprefix $(BINDIR)/, $(BINS))

.PHONY: all clean

all: $(TARGETS)

$(BINDIR):
	mkdir -p $(BINDIR)

# --- 1D core ---

$(BINDIR)/sample_kmers: $(SRCDIR)/sample_kmers.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@ -lz

$(BINDIR)/filter_kmers: $(SRCDIR)/filter_kmers.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@ -lz

$(BINDIR)/compute_structure: $(SRCDIR)/compute_structure.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(EIGEN) $< -o $@

$(BINDIR)/convert_matrix: $(SRCDIR)/convert_matrix.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@ -lz -fopenmp

$(BINDIR)/kgwas: $(SRCDIR)/kgwas.cpp $(SRCDIR)/kgwas_utils.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) -fopenmp -DEIGEN_USE_BLAS $(EIGEN) $< -o $@ -lopenblas -lz

$(BINDIR)/build_index: $(SRCDIR)/build_index.cpp $(SRCDIR)/map_utils.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(VENDOR) -fopenmp $< -o $@

$(BINDIR)/map_to_ref: $(SRCDIR)/map_to_ref.cpp $(SRCDIR)/map_utils.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(VENDOR) -fopenmp $< -o $@

$(BINDIR)/map_merge_filter: $(SRCDIR)/map_merge_filter.cpp $(SRCDIR)/map_utils.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(VENDOR) -fopenmp $< -o $@

# --- 2D epistasis ---

$(BINDIR)/compute_structure_all: $(SRCDIR)/compute_structure_all.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(EIGEN) $< -o $@

$(BINDIR)/kgwas_2d: $(SRCDIR)/kgwas_2d.cpp $(SRCDIR)/kgwas_utils.hpp $(SRCDIR)/chisq.hpp $(SRCDIR)/map_utils.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(VENDOR) -fopenmp -DEIGEN_USE_BLAS $(EIGEN) $< -o $@ -lopenblas -lz

$(BINDIR)/map_2d_pairs: $(SRCDIR)/map_2d_pairs.cpp $(SRCDIR)/map_utils.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(VENDOR) -fopenmp $< -o $@

# --- 2W two-way ---

$(BINDIR)/compute_structure_2w: $(SRCDIR)/compute_structure_2w.cpp $(SRCDIR)/io_2w.hpp $(SRCDIR)/la_2w.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(EIGEN) $< -o $@

$(BINDIR)/compute_projections_2w: $(SRCDIR)/compute_projections_2w.cpp $(SRCDIR)/io_2w.hpp $(SRCDIR)/la_2w.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(EIGEN) $< -o $@

$(BINDIR)/kgwas_2w: $(SRCDIR)/kgwas_2w.cpp $(SRCDIR)/io_2w.hpp $(SRCDIR)/la_2w.hpp $(SRCDIR)/stats_2w.hpp $(SRCDIR)/chisq.hpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) -fopenmp $(EIGEN) $< -o $@

clean:
	rm -rf $(BINDIR)
