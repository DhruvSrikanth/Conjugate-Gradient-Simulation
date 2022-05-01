# Environment variables
PYTHON = python3
GCC = g++-11
MPI = mpiexec
NODES = 4



# Checks
VISUALIZE = "true"
SANITY_CHECK = "false"

# Input to cpp
n = 75

# Run sequence
TARGETS = output sanity_check simulation visualize clean_runtime

# Run the following
target: $(TARGETS)


# Sanity check
sanity_check:
ifeq ($(SANITY_CHECK), "true")
	@echo "Sanity Check:\n"
	$(PYTHON) ./conjugate_gradient.py
	@echo "\n"
endif

output:
	mkdir output

# Compile the following
conjugate_gradient: ./conjugate_gradient.cpp
	$(GCC) ./conjugate_gradient.cpp -o ./conjugate_gradient -I /opt/homebrew/Cellar/open-mpi/4.1.3/include -O3 -ffast-math -march=native -mtune=native -lm -L /opt/homebrew/Cellar/open-mpi/4.1.3/lib -L /opt/homebrew/opt/libevent/lib -lmpi


# Execute the following
simulation: ./conjugate_gradient
	@echo "Simulating...\n"
	$(MPI) --bind-to none -np $(NODES) ./conjugate_gradient $(n)
	@echo "\nSimulation Complete!\n"

# Visualize simulation
visualize:
ifeq ($(VISUALIZE), "true")
	@echo "Visualizing...\n"
	$(PYTHON) ./visualize.py
	@echo "\nVisualization Complete!\n"	
endif

# Clean the following
clean_runtime:
	@echo "Cleaning runtime files...\n"
	rm -rf __pycache__
	rm *.o conjugate_gradient
	rm .DS_Store

clean_text:
	@echo "Cleaning text files...\n"
	rm ./output/*.txt

clean_output:
	@echo "Cleaning output files...\n"
	rm -rf ./output