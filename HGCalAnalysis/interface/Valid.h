#ifndef Validation_HGCalBValidation_Valid_h
#define Validation_HGCalBValidation_Valid_h

#include "TObject.h"
#include <iostream>
#include <vector>

class Valid {
public:
  Valid() : energy_(0), time_(0), timeError_(0), id_(0), flagBits_(0), son_(0) {}
  Valid(float energy, float time, float timeError, uint32_t id, uint32_t flagBits, uint8_t son)
      : energy_(energy), time_(time), timeError_(timeError), id_(id), flagBits_(flagBits), son_(son) {}
  Valid(const Valid &other) {
    energy_ = other.energy_;
    time_ = other.time_;
    timeError_ = other.timeError_;
    id_ = other.id_;
    flagBits_ = other.flagBits_;
    son_ = other.son_;
  }

  virtual ~Valid() {}

  float energy() { return energy_; }
  float time() { return time_; }
  float timeError() { return timeError_; }
  uint32_t id() { return id_; }
  uint32_t flagBits() { return flagBits_; }
  uint8_t son() { return son_; }

  float energy_;       //calibrated energy of the rechit
  float time_;         //time jitter of the UncalibRecHit
  float timeError_;    //time resolution
  uint32_t id_;        //rechit detId
  uint32_t flagBits_;  //rechit flags describing its status (DataFormats/HGCRecHit/interface/HGCRecHit.h)
  uint8_t son_;        //signal over noise

  ClassDef(Valid, 1)
};

typedef std::vector<Valid> ValidCollection;

#endif  //Validation_HGCalBValidation_Valid_h
