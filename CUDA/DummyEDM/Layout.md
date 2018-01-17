# Layout - example of the package layout for CUDA enabled plugins

## Overall __initial__ strategy 
- Complete separation/isolation of CUDA specific code.
- Compiling CUDA code into a shared library (not as an __edm plugin__, although it is also a shared library strictly speaking) together with the rest.
- Kernel invocations is the only thing that requires wrappers

## EDM Plugins (plugins/)
- All the plugins are compiled with g++
- No CUDA specific code (keywords nor kernel invocations) should reside in here

## interface/
- CUDA wrapper __declarations__ should reside in here

## src/
- CUDA wrapper __definitions__ should reside in here
- CUDA Kernels should reside in here
