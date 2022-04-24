# Environment variables
PYTHON = python3
GCC = g++-11

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
conjugate_gradient.o: ./conjugate_gradient.cpp
	$(GCC) -c ./conjugate_gradient.cpp -O3 -ffast-math -march=native -mtune=native -lm

# Link the following
conjugate_gradient: ./conjugate_gradient.o
	$(GCC) ./conjugate_gradient.o -o ./conjugate_gradient

# Execute the following
simulation: ./conjugate_gradient
	@echo "Simulating...\n"
	./conjugate_gradient $(n)
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