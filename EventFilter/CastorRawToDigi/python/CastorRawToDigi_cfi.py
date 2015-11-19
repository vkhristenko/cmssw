import FWCore.ParameterSet.Config as cms

# This version is intended for unpacking standard production data
castorDigis = cms.EDProducer("CastorRawToDigi",
    # Optional filter to remove any digi with "data valid" off, "error" on, 
    # or capids not rotating
    FilterDataQuality = cms.bool(True),
    # Number of the first CASTOR FED.  If this is not specified, the
    # default from FEDNumbering is used.
    CastorFirstFED = cms.int32(690),
    ZDCFirstFED = cms.int32(693),                         
    # FED numbers to unpack.  If this is not specified, all FEDs from
    # FEDNumbering will be unpacked.
<<<<<<< HEAD
    FEDs = cms.untracked.vint32( 690, 691, 692, 693, 722),
=======
    FEDs = cms.untracked.vint32( 690, 691, 692, 693),
>>>>>>> 4e8af02611ac29bc1f616fa2ea121036378a4931
    # Do not complain about missing FEDs
    ExceptionEmptyData = cms.untracked.bool(False),
    # Do not complain about missing FEDs
    ComplainEmptyData = cms.untracked.bool(False),
    # At most ten samples can be put into a digi, if there are more
    # than ten, firstSample and lastSample select which samples
    # will be copied to the digi
    firstSample = cms.int32(0),
    lastSample = cms.int32(9),
    # castor technical trigger processor
    UnpackTTP = cms.bool(True),
    # report errors
    silent = cms.untracked.bool(False),
    #
    InputLabel = cms.InputTag("rawDataCollector"),
    CastorCtdc = cms.bool(False),
    UseNominalOrbitMessageTime = cms.bool(True),
    ExpectedOrbitMessageTime = cms.int32(-1),
    UnpackZDC = cms.bool(True),
<<<<<<< HEAD
=======
    #myMap = cms.VPSet(
     #  cms.PSet(
      ##     object = cms.string('ElectronicsMap'),
        #   file = cms.FileInPath('/afs/cern.ch/user/j/jgomez2/ZDCTest/QuickFix/CMSSW_7_5_0/src/EventFilter/CastorRawToDigi/data/old_ZDC_emap.txt')

     #  ))
>>>>>>> 4e8af02611ac29bc1f616fa2ea121036378a4931

   )
