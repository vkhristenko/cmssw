# Useful Profiling Instructions and Information

## CLI nvprof
- Generate CLI profile `nvprof --print-gpu-trace cmsRun TestGPU/Dummy/python/test_stream_producer_gpu.py`. 

## Generate Metrics for Visual Profiler
- Generate the output file with metrics to be imported into Visual Profiler: `nvprof -o test_newalloc_11122017_0.nvvp cmsRun TestGPU/Dummy/python/test_stream_producer_gpu.py`
