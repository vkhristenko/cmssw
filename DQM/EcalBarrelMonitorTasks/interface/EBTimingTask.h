#ifndef EBTimingTask_H
#define EBTimingTask_H

/*
 * \file EBTimingTask.h
 *
 * $Date: 2007/03/20 12:37:26 $
 * $Revision: 1.3 $
 * \author G. Della Ricca
 *
*/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMServices/Core/interface/MonitorElement.h"

class EBTimingTask: public edm::EDAnalyzer{

public:

/// Constructor
EBTimingTask(const edm::ParameterSet& ps);

/// Destructor
virtual ~EBTimingTask();

protected:

/// Analyze
void analyze(const edm::Event& e, const edm::EventSetup& c);

/// BeginJob
void beginJob(const edm::EventSetup& c);

/// EndJob
void endJob(void);

/// Setup
void setup(void);

/// Cleanup
void cleanup(void);

private:

int ievt_;

edm::InputTag EcalUncalibratedRecHitCollection_;

MonitorElement* meTimeMap_[36];

bool init_;

};

#endif
