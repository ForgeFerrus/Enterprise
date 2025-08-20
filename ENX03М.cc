// ENX03M.cc
// файли Geant4 
#include "G4UserRunAction.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
// файли для роботи з кольорами та візуалізацією
#include "G4UIExecutive.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UImessenger.hh"
#include "G4UImanager.hh"
#include "G4VisAttributes.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4AnalysisManager.hh"
#include "G4GeometryManager.hh"
// файли, що містять класи для геометрії та матеріалів
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
// файли, що містять класи для фізики та частинок
#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4IonConstructor.hh"
// файли для роботи з гістограмами та аналізом
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

namespace ENX03M {
constexpr double target_thickness = 2.0 * mm;
constexpr int NUM_BINS = 512;
constexpr double maxEnergy = 10.0;  // MeV

class SpectrumSD : public G4VSensitiveDetector {
public:
  SpectrumSD(const G4String& name, const G4String& id)
    : G4VSensitiveDetector(name), fId(id), gammaHist(NUM_BINS,0) {}

  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
    auto pt = step->GetPreStepPoint();
    if(pt->GetStepStatus() != fGeomBoundary) return false;
    double E = pt->GetKineticEnergy() / MeV;
    int bin = std::clamp(int(E/maxEnergy*NUM_BINS), 0, NUM_BINS-1);
    if(step->GetTrack()->GetParticleDefinition()->GetParticleName() == "gamma")
      gammaHist[bin]++;
    return true;
  }

  void Save() {
    std::ofstream out(fId + "_gamma.txt");
    for(int i=0; i<NUM_BINS; ++i){
      double E = (i+0.5)*maxEnergy/NUM_BINS;
      out << E << " " << gammaHist[i] << "\n";
    }
  }

private:
  G4String fId;
  std::vector<int> gammaHist;
};

G4LogicalVolume* MakeDetector(G4LogicalVolume* parent, const G4ThreeVector& pos, const G4String& id) {
  auto* vac = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  auto* cyl = new G4Tubs("Tub"+id, 0, 5*cm, 0.5*mm, 0,360*deg);
  auto* lv  = new G4LogicalVolume(cyl, vac, "LV"+id);
  new G4PVPlacement(nullptr, pos, lv, "PV"+id, parent, false,0,true);

  auto* sd = new SpectrumSD("Spec"+id, id);
  G4SDManager::GetSDMpointer()->AddNewDetector(sd);
  lv->SetSensitiveDetector(sd);

  auto* wa = new G4VisAttributes(G4Colour(1,1,1));
  wa->SetForceWireframe(true);
  wa->SetForceSolid(false);
  lv->SetVisAttributes(wa);

  return lv;
}

class Detector_Tungsten : public G4VUserDetectorConstruction {
public:
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
    auto* tungsten = nist->FindOrBuildMaterial("G4_W");

    auto* worldSolid = new G4Box("WorldSolid", 15*cm,15*cm,15*cm);
    auto* worldLog = new G4LogicalVolume(worldSolid, vacuum, "WorldLV");
    auto* worldPhys = new G4PVPlacement(nullptr, {}, worldLog, "WorldPV", nullptr, false, 0);

    auto* targetSolid = new G4Box("TargetSolid", 93/2*mm, 55/2*mm, target_thickness/2);
    auto* targetLog = new G4LogicalVolume(targetSolid, tungsten, "TargetLV");
    new G4PVPlacement(nullptr, {}, targetLog, "TargetPV", worldLog, false, 0);
    targetLog->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));

    MakeDetector(worldLog, G4ThreeVector(0, 0, target_thickness/2 + 1*mm), "D1");

    auto* cylSolid = new G4Tubs("EmptyTube", 0, 5 * cm, 7 * cm, 0, 360 * deg);
    auto* cylLog = new G4LogicalVolume(cylSolid, vacuum, "EmptyTubeLog");
    cylLog->SetVisAttributes(new G4VisAttributes(G4Colour::Grey()));
    G4double cylZ = target_thickness/2 + 1*cm + 7*cm;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, cylZ), cylLog, "EmptyTubePhys", worldLog, false, 0);

    MakeDetector(worldLog, G4ThreeVector(0, 0, target_thickness/2 + 1*cm), "Entry");
    MakeDetector(worldLog, G4ThreeVector(0, 0, target_thickness/2 + 1*cm + 14*cm), "Exit");

    return worldPhys;
  }
};

class PhysicsList : public G4VModularPhysicsList {
public:
  PhysicsList() {
    RegisterPhysics(new G4EmStandardPhysics());
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
  }
  void ConstructParticle() override {
    G4LeptonConstructor().ConstructParticle();
    G4MesonConstructor().ConstructParticle();
    G4BaryonConstructor().ConstructParticle();
    G4IonConstructor().ConstructParticle();
  }
  void SetCuts() override {
    SetCutValue(0.01*mm,"gamma");
    SetCutValue(0.01*mm,"e-");
  }
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun* gun;
public:
  PrimaryGeneratorAction() {
    gun = new G4ParticleGun(1);
    gun->SetParticleMomentumDirection({ 0,0,1 });
    gun->SetParticlePosition({ 0,0,-10 * cm });
    gun->SetParticleEnergy(8.7 * MeV);
    gun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("e-"));
  }
  ~PrimaryGeneratorAction() override { delete gun; }
  void GeneratePrimaries(G4Event* evt) override {
    gun->GeneratePrimaryVertex(evt);
  }
};

class RunAction : public G4UserRunAction {
public:
  void EndOfRunAction(const G4Run*) override {
    auto* SDM = G4SDManager::GetSDMpointer();
    for (const auto& id : { "D1", "Entry", "Exit" }) {
      auto* sd = dynamic_cast<SpectrumSD*>(SDM->FindSensitiveDetector("Spec"+id));
      if(sd) sd->Save();
    }
    G4cout << "[Run] Spectra saved.\n";
  }
};

class ActionInitialization : public G4VUserActionInitialization {
public:
  void Build() const override {
    SetUserAction(new PrimaryGeneratorAction());
    SetUserAction(new RunAction());
  }
};

} // namespace ENX03M
// ─────────────────────────────────────────────────────
  // Головна функція запуску задачі
// ===== Інтеграція з Qt-мостом або main() =====
extern "C" void InitUserClasses(G4RunManager* run) {
  run->SetUserInitialization(new ENX03M::Detector_Tungsten());
  run->SetUserInitialization(new ENX03M::PhysicsList());
  run->SetUserInitialization(new ENX03M::ActionInitialization());
}
namespace Enterprise {
    int StartFromBridge(int argc, char** argv);  // або без namespace, якщо нема
}
int main(int argc, char** argv) {
  return Enterprise::StartFromBridge(argc, argv);
}



