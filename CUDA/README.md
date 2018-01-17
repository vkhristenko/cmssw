# CMSSW CUDA package

## Set up
- Get a __DEVEL__ cmssw version
```
cmsrel CMSSW_10_0_DEVEL_X_2017-11-02-2300
cd CMSSW_10_0_DEVEL_X_2017-11-02-2300
```
- `scram setup cuda`
- `cmsenv`
- Pull the branch with cuda samples and edm dummy examples
```
git cms-merge-topic vkhristenko:testgpu
scram -b -v -j 8
```

## Samples
- Original CUDA Nvidia Samples
- Executables reside in `$CMSSW_BASE/test/$SCRAM_ARCH/`

## DummyEDM
- Example of using GPU offload with CMSSW EDM Stream Producers
- For further info see DummyEDM/README.md or DummyEDM/Layout.md
